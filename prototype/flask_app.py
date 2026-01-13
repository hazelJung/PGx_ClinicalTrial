"""
Flask Web Application - Virtual Population PBPK Simulator
========================================================
HTML/CSS 파일을 직접 수정할 수 있는 Flask 기반 웹앱

구조:
- templates/index.html : 메인 HTML 파일
- static/css/style.css : 스타일시트
- static/js/app.js : JavaScript
"""

from flask import Flask, render_template, jsonify, request
import numpy as np
import requests

# 같은 폴더의 모듈 import
from engine import PopulationGenerator, Ethnicity
from pbpk_model import (
    DrugParameters, SimulationConfig, PhysiologicalParameters,
    run_population_simulation
)


def fetch_pubchem_data(drug_name: str) -> dict:
    """PubChem API에서 약물 정보 가져오기"""
    try:
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        props_url = f"{base_url}/compound/name/{drug_name}/property/MolecularWeight,XLogP,IUPACName/JSON"
        
        response = requests.get(props_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            
            return {
                'found': True,
                'name': drug_name,
                'mw': props.get('MolecularWeight', 0),
                'log_p': props.get('XLogP', 2.0),
                'iupac_name': props.get('IUPACName', '')
            }
        else:
            return {'found': False, 'name': drug_name}
            
    except Exception as e:
        return {'found': False, 'name': drug_name, 'error': str(e)}


app = Flask(__name__)


# ============================================================================
# 페이지 라우트
# ============================================================================

@app.route('/')
def index():
    """메인 페이지"""
    return render_template('index.html')


# ============================================================================
# API 엔드포인트
# ============================================================================

@app.route('/api/generate-population', methods=['POST'])
def api_generate_population():
    """Step 1: 가상 인구집단 생성 API"""
    try:
        data = request.json
        
        # 민족 비율 정규화
        total_eth = data['eth_asian'] + data['eth_european'] + data['eth_african']
        if total_eth == 0:
            total_eth = 100
        
        ethnicity_dist = {
            Ethnicity.EAST_ASIAN: data['eth_asian'] / total_eth,
            Ethnicity.EUROPEAN: data['eth_european'] / total_eth,
            Ethnicity.AFRICAN: data['eth_african'] / total_eth
        }
        
        # 인구집단 생성
        generator = PopulationGenerator(
            n_subjects=data['n_subjects'],
            ethnicity_distribution=ethnicity_dist,
            age_range=(data['age_min'], data['age_max']),
            gender_ratio=data['gender_ratio'] / 100,
            weight_mean=data['weight_mean'],
            weight_sd=data['weight_sd'],
            base_cl_int=data['base_clint'],
            random_seed=None  # 매번 새로운 인구집단
        )
        
        population = generator.generate()
        pop_summary = generator.get_population_summary(population)
        
        # 개인별 데이터 (프론트엔드용)
        individuals = []
        metabolizer_map = {
            'Poor Metabolizer (PM)': 'PM',
            'Intermediate Metabolizer (IM)': 'IM',
            'Normal Metabolizer (NM)': 'NM',
            'Ultra-rapid Metabolizer (UM)': 'UM'
        }
        
        for ind in population:
            individuals.append({
                'id': ind.subject_id,
                'age': ind.age,
                'gender': ind.gender,
                'weight': round(ind.weight, 1),
                'height': round(ind.height, 1),
                'bmi': round(ind.bmi, 1),
                'ethnicity': ind.ethnicity.value,
                'metabolizer': metabolizer_map.get(ind.metabolizer_status.value, 'NM'),
                'activity_score': round(ind.combined_activity_score, 2),
                'cyp2c19': '/'.join(ind.cyp2c19_genotype),
                'cyp3a4': '/'.join(ind.cyp3a4_genotype),
                'phys_params': {
                    'body_weight': ind.phys_params.body_weight,
                    'v_plasma': ind.phys_params.v_plasma,
                    'v_liver': ind.phys_params.v_liver,
                    'q_liver': ind.phys_params.q_liver,
                    'cl_int': ind.phys_params.cl_int,
                    'cl_renal': ind.phys_params.cl_renal,
                    'activity_score': ind.phys_params.activity_score
                }
            })
        
        return jsonify({
            'success': True,
            'individuals': individuals,
            'summary': pop_summary
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/api/run-simulation', methods=['POST'])
def api_run_simulation():
    """Step 2: 약물 시뮬레이션 API"""
    try:
        data = request.json
        
        # 약물 파라미터
        drug_params = DrugParameters(
            name=data.get('drug_name', 'Unknown'),
            log_p=data['log_p'],
            f_u=data['f_u'],
            v_d=data['v_d'],
            k_a=data['k_a'],
            f=data['bioavail']
        )
        
        # 시뮬레이션 설정
        sim_config = SimulationConfig(
            dose=data['dose'],
            route='oral',
            t_max=24,
            n_points=241
        )
        
        # 인구집단 생리학적 파라미터 복원
        population_phys = []
        for ind in data['population']:
            phys = ind.get('phys_params', {})
            population_phys.append(PhysiologicalParameters(
                body_weight=phys.get('body_weight', 70),
                v_plasma=phys.get('v_plasma', 3.0),
                v_liver=phys.get('v_liver', 1.5),
                q_liver=phys.get('q_liver', 90),
                cl_int=phys.get('cl_int', 10),
                cl_renal=phys.get('cl_renal', 0),
                activity_score=phys.get('activity_score', 1.0)
            ))
        
        # 시뮬레이션 실행
        results = run_population_simulation(drug_params, population_phys, sim_config)
        
        return jsonify({
            'success': True,
            'time': results['time'].tolist(),
            'mean_concentration': results['mean_concentration'].tolist(),
            'ci_lower': results['ci_lower'].tolist(),
            'ci_upper': results['ci_upper'].tolist(),
            'individual_curves': results['individual_curves'][:50].tolist(),
            'cmax_distribution': results['cmax_distribution'].tolist(),
            'auc_distribution': results['auc_distribution'].tolist()
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/api/fetch-pubchem', methods=['GET'])
def api_fetch_pubchem():
    """PubChem에서 약물 정보 가져오기"""
    drug_name = request.args.get('drug_name', '')
    
    if not drug_name:
        return jsonify({'found': False, 'error': '약물명을 입력하세요'})
    
    data = fetch_pubchem_data(drug_name)
    return jsonify(data)


# ============================================================================
# 서버 실행
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*60)
    print("  🧬 Virtual Population PBPK Simulator (Flask Version)")
    print("="*60)
    print("  서버 시작 중...")
    print("  브라우저에서 http://localhost:5000 을 열어주세요")
    print("="*60)
    print("\n  📁 파일 구조:")
    print("     templates/index.html  - HTML 수정")
    print("     static/css/style.css  - CSS 수정")
    print("     static/js/app.js      - JavaScript 수정")
    print("="*60 + "\n")
    
    app.run(debug=True, port=5000)
