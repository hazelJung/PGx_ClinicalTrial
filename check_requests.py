import requests
import time

print("1. requests 모듈 임포트 시작...")
start_time = time.time()
try:
    import requests
    print(f"   - 임포트 성공! (소요시간: {time.time() - start_time:.4f}초)")
except ImportError:
    print("   - 임포트 실패! 가상환경 설정을 확인해주세요.")
    exit()

print("\n2. PubChem API 연결 테스트 중...")
cid = "11962412"
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"

try:
    # 5초 타임아웃 설정
    response = requests.get(url, timeout=5)
    print(f"   - 응답 코드: {response.status_code}")
    if response.status_code == 200:
        print("   - 데이터 로드 성공!")
        data = response.json()
        print(f"   - 데이터 크기: {len(str(data))} bytes")
    else:
        print("   - 데이터 로드 실패")
except requests.exceptions.Timeout:
    print("   - [에러] 요청 시간이 초과되었습니다. (Timeout)")
except Exception as e:
    print(f"   - [에러] 연결 중 문제 발생: {e}")
