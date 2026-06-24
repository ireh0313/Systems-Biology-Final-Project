# Fig. 4/5 최종 재현 요약

## 최종 설정

최종 Fig. 4/5 재현은 corrected 4TF_self 모델로 진행했다.

```text
TF: ATF1, ATF7, CREB1, CREM
Competition: false
SelfCompetition: true
full fitting init_loop: 100000
real seed 반복: 5회
scramble control: 10회
```

여기서 중요한 점은 transcpp의 `Competition` 옵션이 논문에서 말하는 TFBS competitive binding/self-competition과 같은 개념이 아니라는 것이다. 그래서 최종 분석에서는 `Competition=false`로 두고, 4TF_self 모델에 해당하는 `SelfCompetition=true`를 사용했다.

## 대표 Seed

최종 그림에는 `seed0`을 대표 모델로 사용했다.

```text
seed0 single-hit r = 0.834
seed0 Fig. 4J multi-hit r = 0.764
seed0 Fig. 5 A/T profile r = 0.869
```

seed0은 Fig. 4J multi-hit validation 성능이 5개 seed 중 가장 좋고, single-hit 및 Fig. 5 성능도 높은 편이라 대표 그림으로 적합하다.

## Fig. 4에서 보여줄 내용

Fig. 4 최종 plate는 다음을 보여준다.

```text
A: single-hit mutation profile, observed vs predicted
B: single-hit observed-predicted scatter
C: multi-hit validation scatter, Fig. 4J-style
D: real seed 5회와 scramble 10회 비교
```

핵심 해석은 real label로 학습한 모델이 scramble control보다 훨씬 높은 correlation을 보인다는 것이다.

```text
real seed single-hit mean r = 0.828
scramble true-label mean r = -0.093
scramble true-label max r = 0.280
```

즉 모델이 무작위 label을 단순히 fitting한 것이 아니라, 실제 sequence와 enhancer activity 사이의 관계를 학습했다고 해석할 수 있다.

## Fig. 5에서 보여줄 내용

Fig. 5 최종 plate는 다음을 보여준다.

```text
A: A/T substitution activity profile
B: CRE1/CRE4 주변 TFBS occupancy map
C: fitted activation coefficients
D: CRE1/CRE4 contribution proxy
```

Fig. 5 A/T substitution profile은 real seed 5회 평균으로도 안정적이었다.

```text
Fig. 5 A/T profile mean r = 0.869
range = 0.863-0.876
```

CRE1과 CRE4는 mutation profile에서 activity 감소가 크게 나타나는 주요 구간이고, occupancy/contribution 분석은 이 구간의 TF binding이 모델 예측에 중요하다는 점을 설명하기 위한 시각화이다.

## 최종 제출 파일

가장 중요한 파일은 아래 두 개이다.

```text
final_fig4_competition_off.svg
final_fig5_competition_off.svg
```

보고서에 수치나 method를 넣을 때는 아래 파일을 보면 된다.

```text
fig45_final_metrics.csv
fig45_final_submission_README.md
fig45_corrected_methods_for_report.md
```

전체 반복 결과 원본은 아래 폴더에 있다.

```text
outputs/fig45_competition_off_repeats/
```
