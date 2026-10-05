# Exomiser와 LIRICAL 및 AVI 분석 코드

제출 이름은 **hideonbocchi**다. 이 저장소는 v1.1.0 코드를 실행 단계별로 정리한 **v1.1.1**이며, 후보·순위·Tier·EPCR을 유지한다. 제출용 CSV에는 후보별 notes만 추가했다.

- `scripts/`: 실행 스크립트. 정규화, 주석, read/phase, AVI는 하위 폴더로 구분했다.
- `src/`: 공통 분석 로직. `canonical/`의 기존 판정 함수는 원본을 유지했다.
- `config/`, `examples/`: 공개 설정과 입력 경로 예시.
- `tests/`: 합성 자료와 공개 대조군을 이용한 검증.
- `docs/`: 방법, 입력 형식, 분석 단계, 한계.
- `provenance/`: 원본 출처, 파일 이동 기록, 검증 결과.

## 재현 범위

검증된 재현은 저장된 Exomiser·LIRICAL 출력과 주석·read·검토 관측에서 시작해, 120쌍의 판정과 제출 CSV를 다시 만든다. AVI는 저장 응답을 읽으므로 재조회나 API 키가 필요 없다. 원본 WGS부터 모든 대규모 분석을 다시 실행한 것은 아니다. 문헌 검토는 별도의 기록된 입력이다.

## 실행

Python 3.11, bcftools 1.22, 접근 권한이 있는 비공개 입력 묶음과 기준 유전체가 필요하다. 저장소 루트에서 실행하고 `/private/...`를 실제 경로로 바꾼다. 출력 폴더는 새 경로여야 한다.

```bash
export TMPDIR=/private/scratch
python -B scripts/reproduce.py --config /private/reproduction.private.json --avi-config /private/avi.private.json --output /private/new-replay
python -B scripts/verify_reproduction.py --result /private/new-replay --reference /private/verification_only
```

AVI를 생략하려면 `--avi-config` 옵션만 빼면 된다. 기본 검증과 AVI 검증의 PASS를 함께 확인한다. 전체 분석 단계는 [UPSTREAM](docs/UPSTREAM.md), 세부 사용법은 [영문 README](README.md)를 따른다.

## 제출

포털에는 비공개 CSV `hideonbocchi_pairprior_v1.csv`와 [영문 보고서](METHODS.md)를 `hideonbocchi_track1_report.md`로 저장해 사용한다. Anthropic Claude(Team plan for scientists)는 분석을 지원했고, OpenAI Codex(ChatGPT Pro 20x)는 코드 작성과 보고서 작성을 지원했다. 제출 자료는 모델 학습에 사용하지 않았다(참가자 확인). 데이터 보관 설정은 미확인이다. 환자 자료, CSV, API 응답·키는 GitHub에 올리지 않는다. 이 저장소 게시와 포털 제출은 별개다.

제출용 CSV는 후보별 근거와 한계를 notes에 추가한 검토본이다. 다른 열과 행 순서는 원본 재현 출력과 같다. 재현 코드와 바이트 비교 검증기는 기존 `Phase unresolved.` CSV를 계속 생성·검증한다. notes 보완은 선정 규칙이 아닌 제출 문구 편집 단계이며, 포털에는 검토본 CSV를 사용한다.
