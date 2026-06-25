from __future__ import annotations

import math
import re
import statistics
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

WT_ID = "synCRE_Promega_0"
MUTATION_RE = re.compile(r"^scanmut_single_pos_(\d+)_([ACGT])$")


def _column_index(cell_reference: str) -> int:
    letters = "".join(char for char in cell_reference if char.isalpha())
    value = 0
    for char in letters.upper():
        value = value * 26 + ord(char) - ord("A") + 1
    return value - 1


def read_xlsx_sheet(path: Path, sheet_name: str) -> list[list[object]]:
    spreadsheet_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    document_rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    package_rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"

    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{spreadsheet_ns}}}si"):
                shared_strings.append(
                    "".join(node.text or "" for node in item.iter(f"{{{spreadsheet_ns}}}t"))
                )

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rel_id = None
        for sheet in workbook.iter(f"{{{spreadsheet_ns}}}sheet"):
            if sheet.attrib.get("name") == sheet_name:
                rel_id = sheet.attrib.get(f"{{{document_rel_ns}}}id")
                break
        if rel_id is None:
            raise ValueError(f"Sheet {sheet_name!r} was not found in {path}")

        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = None
        for rel in rels.iter(f"{{{package_rel_ns}}}Relationship"):
            if rel.attrib.get("Id") == rel_id:
                target = rel.attrib.get("Target")
                break
        if target is None:
            raise ValueError(f"Could not resolve sheet {sheet_name!r} in {path}")

        sheet_path = target.lstrip("/")
        if not sheet_path.startswith("xl/"):
            sheet_path = f"xl/{sheet_path}"
        sheet_root = ET.fromstring(archive.read(sheet_path))

    rows: list[list[object]] = []
    for row in sheet_root.iter(f"{{{spreadsheet_ns}}}row"):
        parsed: dict[int, object] = {}
        for cell in row.findall(f"{{{spreadsheet_ns}}}c"):
            idx = _column_index(cell.attrib.get("r", "A1"))
            value_node = cell.find(f"{{{spreadsheet_ns}}}v")
            cell_type = cell.attrib.get("t")
            if cell_type == "inlineStr":
                value: object = "".join(
                    node.text or "" for node in cell.iter(f"{{{spreadsheet_ns}}}t")
                )
            elif value_node is None:
                value = ""
            elif cell_type == "s":
                value = shared_strings[int(value_node.text or 0)]
            else:
                text = value_node.text or ""
                try:
                    value = float(text)
                except ValueError:
                    value = text
            parsed[idx] = value
        if parsed:
            rows.append([parsed.get(i, "") for i in range(max(parsed) + 1)])
    return rows


def read_measured_activity(path: Path) -> dict[str, float]:
    rows = read_xlsx_sheet(path, "single-hit")
    headers = [str(value).strip().lower() for value in rows[0]]
    name_idx = headers.index("name")
    activity_idx = headers.index("expression level")
    measured: dict[str, float] = {}
    for row in rows[1:]:
        if len(row) <= max(name_idx, activity_idx):
            continue
        variant_id = str(row[name_idx]).strip()
        if variant_id:
            measured[variant_id] = float(row[activity_idx])
    if WT_ID not in measured:
        raise ValueError(f"WT row {WT_ID!r} was not found in {path}")
    return measured


def read_rates(path: Path) -> dict[str, float]:
    rates: dict[str, float] = {}
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].lower() == "id":
                continue
            value = float(fields[-1])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Invalid rate in {path}: {value}")
            rates[fields[0]] = value
    if WT_ID not in rates:
        raise ValueError(f"WT row {WT_ID!r} was not found in {path}")
    return rates


def pearson_r(x_values: list[float], y_values: list[float]) -> float:
    x_mean = statistics.fmean(x_values)
    y_mean = statistics.fmean(y_values)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    x_ss = sum((x - x_mean) ** 2 for x in x_values)
    y_ss = sum((y - y_mean) ** 2 for y in y_values)
    return numerator / math.sqrt(x_ss * y_ss)


def model_metrics(rates: dict[str, float], measured: dict[str, float]) -> tuple[float, float, int]:
    measured_wt = measured[WT_ID]
    predicted_wt = rates[WT_ID]
    observed: list[float] = []
    predicted: list[float] = []
    for variant_id, value in measured.items():
        if MUTATION_RE.match(variant_id) and variant_id in rates:
            observed.append(math.log2(value / measured_wt))
            predicted.append(math.log2(rates[variant_id] / predicted_wt))
    r_value = pearson_r(observed, predicted)
    rmse = math.sqrt(
        statistics.fmean((pred - obs) ** 2 for obs, pred in zip(observed, predicted))
    )
    return r_value, rmse, len(observed)
