# 02. Kim2013 example fitting check

This step verifies that the compiled `transcpp` executable can read and fit the bundled example XML.

The example file is included in the original `transcpp` repository:

```text
~/transcpp-master/fits/kim2013.xml
```

Run:

```bash
bash run_kim2013_smoke_test.sh
```

Expected success signal:

```text
Created Nuclei
output verified!
<Output>
```

The Boost `bind.hpp` deprecation note can appear while building `unfold`; it is a warning, not a failure.

