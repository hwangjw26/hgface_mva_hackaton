# Exomiser–LIRICAL 우선순위 분석과 AVI 보조 주석

Participant: **hideonbocchi** · [GitHub repository](https://github.com/hwangjw26/hgface_mva_hackaton)
Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1

## 접근법과 자료

제공된 단일 환자 WGS VCF와 환자 HPO 7개를 분석하고 가족력은 별도로 검토했다. 주최 측은 정답의 compound-heterozygous 구조를 공개했다. 증례 자료는 접근 제한 대상이며 보조 데이터베이스와 문헌은 공개 자료다.

Exomiser 15.1.0과 hg38/phenotype 2602의 coding/splice·genome preset에서 두 이형접합 allele이 기여하는 native AR 유전형 120쌍을 확보했다. 설정별 원 순위를 유지했고 알려진 MVA 유전자로 탐색을 제한하지 않았다.

LIRICAL 2.4.1로 전체 VCF의 기준 분석과 HPO 하나씩 제외한 7회 분석을 수행했다. non-global, pathogenicity threshold 0.8, 기본 genotype penalty를 사용했고 유전자·질환 제한과 `--sdwndv`는 없었다. 전체 JSON에서 30쌍의 AR 질환 문맥을 확인했다. 질환 반환이나 native allele 적격성은 해당 쌍의 원인성을 입증하지 않는다.

자동 주석과 출처 기반 검토로 정확 allele의 ClinVar 기록(2026-08-22), transcript, 집단 빈도, 생물학적 문헌, read와 phase를 평가했다. bcftools 1.22, VEP 115/gnomAD v4.1 주석, HPO 2026-06-23, HGNC 2026-04-07을 사용했다. MVA를 포함한 유전자 수준 문헌과 allele 수준 지지를 구분하고 근거 결측은 미확인으로 유지했다.

## AVI 보조 분석

2026-10-05 Google AlphaGenome Atlas SDK 0.9.0으로 SNV 201개의 AVI 점수와 18개 feature 기여도를 받았다. 승인된 조회는 GRCh38 chromosome·position·REF·ALT만 전송했다. Indel 37개는 미평가로 남겼다.

직접 API quantile을 [Qu와 Hakonarson 사전 논문](https://doi.org/10.64898/2026.09.26.754650)의 Table 2 기준 **>=0.99379**와 비교했다. 두 변이 모두 SNV인 93쌍 중 양쪽 기준 충족 4쌍, 한쪽 18쌍, 둘 다 미달 71쌍이며 27쌍은 불완전 비교였다. 전체 120쌍에 AVI를 주석으로 연결했고 제외·재순위·제출 값 변경에는 사용하지 않았다. AlphaMissense·보존도 등 공유 입력과 ClinVar 기반 평가 때문에 독립 근거로 가정하지 않았다. Quantile은 병원성 확률이 아니며 feature 기여도는 기능 실험이 아니다.

## 선정과 한계

잠정 기전 검토 tier를 충족한 4쌍을 선정했고 더 강한 allele 근거 tier를 충족한 쌍은 없었다. 두 Exomiser 설정 중 더 높은 AR 순위, exome 순위, genome 순위, 식별자 순으로 정렬했다. 네 행은 주된 질환의 대안 가설인 `primary`이며 부수적 소견은 제출하지 않는다. 비공개 정답이나 leaderboard 피드백을 사용하지 않았다.

EPCR 1.00·0.75·0.50·0.25는 순위 표현값이며 보정된 원인 확률이 아니다. 네 쌍 모두 phase 미정이며 3쌍은 양 allele의 강한 read 지지가 있다. 연구 가설이며 구조·반복·모자이크 기전은 포괄적으로 평가하지 않았다.

## 재현성과 자원

v1.1.1은 기존 v1.1.0 코드의 폴더 구조를 정리하고 판단을 유지한다. Native 출력과 명시된 주석·검토 관측에서 기준 분석을 재구성한 뒤 선택적으로 hash를 고정한 응답에서 AVI 주석을 추가하고 독립 검증한다. 저장 응답 재현에는 API 키나 네트워크 요청이 필요 없다. 서비스는 데이터 snapshot 버전을 제공하지 않는다. 패키징을 위해 대규모 앞 단계 분석은 반복하지 않았다.

관측 wall time은 Exomiser 5분47초/10분28초, LIRICAL 8회 34분32초, 표적 phase 확장 86분38초다. 다른 작업은 포함하지 않으며 총 노력과 비용은 기록하지 않았다.

**AI 지원:** Anthropic Claude(Team plan for scientists)는 분석을 지원했고, OpenAI Codex(ChatGPT Pro 20x)는 코드 작성과 보고서 작성을 지원했다. 계정의 학습·보관 설정은 미확인이다.

제출 CSV의 notes에는 검토한 후보별 근거를 추가했으며 다른 열은 변경하지 않았다.

자료를 제공한 가족과 대회 주최 측에 감사한다.
