# 7dt-hips — 7DT colour HiPS (RIS DP0 first)

7DT single processed science frame(v2)으로 밴드별 FITS HiPS를 만들고, CDS
Thomas Boch의 색 합성 레시피(기본값 그대로)로 colour PNG HiPS를 만든다.
첫 시험은 **RIS DP0 7개 IMS 타일**, 실행은 **skynet**(`/lyman` 마운트됨,
`/balmer`는 안 됨).

## 단계

| 단계 | 스크립트 | 입력 → 출력 |
|---|---|---|
| 1 | `hips7dt/frame_index.py` | v2 frame 헤더 → `index/dp0_frames.csv` |
| 2 | `hips7dt/select_link.py` | index → `hips_in/dp0/<band>/` symlink 트리 + `manifest.csv`, `zp_report.csv` |
| 3 | `hips7dt/run_hipsgen.sh` | symlink 트리 → `hips/dp0/<band>/` FITS HiPS (Hipsgen) |
| 4 | `hips7dt/make_color_hips.py` | `hips/dp0/m*/` → `hips/dp0/color/` PNG HiPS |

모든 Python 단계에 `--demo`(데이터 없이 자체 검사)가 있다.
`make_color_hips.py --demo`는 Thomas 원본(`docs/thomas/`)과 같은 입력에서
**픽셀 단위로 동일한 타일**이 나오는지 확인한다.

## skynet 준비 (한 번)

```bash
pip install astropy numpy pandas mocpy pillow
mkdir -p ~/bin && wget -O ~/bin/Hipsgen.jar https://aladin.cds.unistra.fr/java/Hipsgen.jar
java -jar ~/bin/Hipsgen.jar -h | less     # 파라미터 이름 확인 (아래 주의)
python3 hips7dt/frame_index.py --demo && python3 hips7dt/select_link.py --demo \
  && python3 hips7dt/make_color_hips.py --demo
```

## DP0 실행 순서

```bash
# 0) v2 frame 하나의 헤더 확인 -> ZP/seeing/depth 키 이름
python3 hips7dt/frame_index.py --dump-header <v2 single frame 하나>

# 1) 인덱스 ({tile}은 타일별로 펼쳐짐, coadd/ 경로는 자동 제외)
python3 hips7dt/frame_index.py --glob '<v2 root>/*/{tile}/*/*.fits' \
    --tiles config/dp0_tiles.txt --out index/dp0_frames.csv --workers 4

# 2) 선택 (타일×밴드마다 seeing 상위 10장) - 먼저 dry-run으로 숫자 확인
python3 hips7dt/select_link.py --index index/dp0_frames.csv --out hips_in/dp0 \
    --tiles config/dp0_tiles.txt --exclude-dates 2025-09-04 \
    --per-tile-band 10 --require-zp --dry-run
#    결과가 괜찮으면 --dry-run 빼고 다시

# 3) 밴드 하나로 Hipsgen 시험 -> 괜찮으면 전 밴드
BANDS=m650 IN=hips_in/dp0 OUT=hips/dp0 bash hips7dt/run_hipsgen.sh
IN=hips_in/dp0 OUT=hips/dp0 bash hips7dt/run_hipsgen.sh

# 4) 색 합성
python3 hips7dt/make_color_hips.py --in-root hips/dp0 --out hips/dp0/color --workers 4
```

확인: Aladin Desktop에서 `hips/dp0/color` 폴더를 열거나, 로컬 웹서버
(`python3 -m http.server -d hips/dp0`)에 Aladin Lite를 붙인다.

## 주의할 점

- **Hipsgen 파라미터**는 이 저장소를 만든 클라우드 세션에서 실행해 보지
  못했다 (`in= out= id= order= bitpix= maxThread= skyval=auto`). `-h` 출력과
  다르면 `run_hipsgen.sh`를 고친다. 색 스크립트는 각 밴드 `properties`의
  `hips_pixel_cut`과 `Norder3/Allsky.fits`, `Moc.fits`를 필요로 한다.
- **symlink**: Hipsgen이 디렉터리 symlink가 아닌 *파일* symlink를 따라가는지
  m650 한 밴드로 먼저 확인한다 (로그에 읽은 파일 수).
- **ZP 산포**: Hipsgen은 파일마다 스케일을 다르게 줄 수 없다. `zp_report.csv`의
  밴드별 `zp_range`가 0.1 mag를 넘으면 (unit별 감도 차이) 이음매가 생긴다.
  그때는 `--scaled-copy 25` (공통 ZP로 스케일한 float32 사본; 용량 출력됨)를
  시험한다. ZP가 노출시간을 포함하는 정의인지도 헤더에서 확인할 것.
- **DP0는 날짜가 많다** (타일당 140–170일). 전부 넣지 말고 `--per-tile-band`로
  자른다. 10장 × 7타일 × 20밴드 ≈ 1,400장.
- **DP0 감사에서 나온 문제 영상**(2025-06/07 m425·m625·m850, 2025-09-04 전 밴드,
  2025-10 적색 밴드)은 coadd 기준 목록이다. v2에서 고쳐졌는지 모르므로 처음엔
  2025-09-04만 빼고, 결과에서 해당 밴드를 확인한다. proton의
  `outputs/release_integrity/{bands_dp0_culprits,epoch_bad_images}.csv`를
  복사해 오면 `--exclude-images a.csv,b.csv`로 뺄 수 있다.
- **용량**: 7타일 ≈ 8–9 deg², order 10 FITS HiPS ≈ 0.43 GB/deg²/밴드 →
  밴드당 ~4 GB, 20밴드 ~80 GB.

## 개선 작업 (기본값으로 결과를 본 다음)

1. gradient: `skyval=auto` 유무 비교, v2 frame이 배경 차감 상태인지 확인
2. ZP 통일: `--scaled-copy` 결과와 비교해 이음매가 줄어드는지
3. 밴드 결측 경계의 색 띠, order별 밝기 튐 (Thomas의 order blackpoint)
4. 같은 필드를 coadd 입력으로 만든 HiPS와 비교 (우주선·위성 궤적)
