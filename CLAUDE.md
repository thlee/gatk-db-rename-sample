# gatk-db-rename-sample

GATK GenomicsDB에서 샘플 이름을 변경하는 도구. DB 재구축 없이 callset.json 메타데이터를 직접 수정.

## 프로젝트 상태

- 초기 릴리즈 완료 (v0.1.0)
- GitHub: https://github.com/thlee/gatk-db-rename-sample
- 라이선스: MIT

## 구조

```
gatk_db_rename_sample.py   # 메인 스크립트 (단일 파일, 표준 라이브러리만 사용)
README.md                  # 영문 사용 설명서
LICENSE                    # MIT
```

## 기능

- `list` — DB 내 샘플 목록/수 조회
- `rename` — 개별(`--old/--new`) 또는 일괄(`--map`) 이름 변경
- `validate` — 모든 interval workspace 간 샘플 일치 검증
- `restore` — 백업에서 복원

## 핵심 원리

- GenomicsDB의 `callset.json` 내 `sample_name` 필드만 수정
- `row_idx`는 절대 변경하지 않음 (TileDB 배열 내 데이터 위치)
- 엔트리 삭제 불가 (TileDB에 row 데이터가 남아있어 GATK 오류 발생)
- 삭제가 필요하면 `__REMOVED__` 접두사로 rename 후 bcftools로 필터링

## 개발 참고

- Python 3.6+ 호환 유지 (f-string 사용, walrus operator 미사용)
- 외부 의존성 없음 (표준 라이브러리만)
- 저자: TH Lee (기획/디렉션/리뷰), Claude/Anthropic (코드)

## 테스트

실제 GenomicsDB에서 테스트 완료:
- rename: callset.json 수정 후 GenotypeGVCFs 정상 동작 확인
- remove(엔트리 삭제): GATK 오류 발생 확인 → rename 방식으로 대체
- 2,639 샘플, 142 interval workspace 환경에서 검증
