"""
Virtual Population PBPK Simulation Platform
=============================================
임상 1상 전환을 위한 가상 인구집단 PBPK 시뮬레이션 플랫폼

Features:
    - 가상 인구집단 생성 (1,000+ 명)
    - 다중 구획 PBPK 모델 시뮬레이션
    - 실시간 PK 곡선 시각화
    - 안전 마진 분석

Author: PBPK Simulation Team
"""

import numpy as np
import requests
from dash import Dash, html, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from pbpk_model import (
    PBPKModel, DrugParameters, PhysiologicalParameters, 
    SimulationConfig, run_population_simulation
)
from engine import (
    PopulationGenerator, Ethnicity, MetabolizerStatus,
    calculate_safety_margin
)

# ============================================================================
# App Initialization
# ============================================================================

app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    title="PBPK Simulator",
    suppress_callback_exceptions=True
)

# Color Palette
COLORS = {
    'primary': '#1E3A5F',      # Deep Medical Blue
    'secondary': '#2ECC71',    # Safety Green
    'warning': '#F39C12',      # Caution Orange
    'danger': '#E74C3C',       # Alert Red
    'background': '#F8F9FA',   # Clean White-Gray
    'card': '#FFFFFF',
    'text': '#2C3E50',
    'text_muted': '#7F8C8D',
    'border': '#E5E8EB',
    'gradient_start': '#667eea',
    'gradient_end': '#764ba2'
}

# ============================================================================
# PubChem API Integration
# ============================================================================

def fetch_pubchem_data(drug_name: str) -> dict:
    """PubChem에서 약물 데이터 가져오기"""
    try:
        # Search for compound
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/MolecularWeight,XLogP,IUPACName/JSON"
        response = requests.get(search_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            return {
                'found': True,
                'name': drug_name,
                'mw': props.get('MolecularWeight', 300),
                'log_p': props.get('XLogP', 2.0),
                'iupac_name': props.get('IUPACName', '')
            }
    except Exception as e:
        print(f"PubChem API error: {e}")
    
    return {'found': False}

# ============================================================================
# UI Components
# ============================================================================

def create_header():
    """헤더 컴포넌트"""
    return dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.Div([
                        html.Span("🧬", style={'fontSize': '2rem', 'marginRight': '10px'}),
                        html.Span("Virtual Population PBPK Simulator", 
                                 style={'fontSize': '1.5rem', 'fontWeight': '600', 'color': 'white'})
                    ], style={'display': 'flex', 'alignItems': 'center'})
                ], width='auto'),
            ], align='center', className='w-100'),
        ], fluid=True),
        color=COLORS['primary'],
        dark=True,
        className='mb-4',
        style={'boxShadow': '0 2px 10px rgba(0,0,0,0.1)'}
    )


def create_population_card():
    """인구집단 설정 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("👥 인구집단 설정", className='mb-0', 
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            # Total N
            html.Label("총 인원수 (N)", className='fw-semibold text-muted mb-1'),
            dcc.Slider(
                id='n-subjects',
                min=100, max=2000, step=100, value=1000,
                marks={100: '100', 500: '500', 1000: '1K', 1500: '1.5K', 2000: '2K'},
                tooltip={'placement': 'bottom', 'always_visible': True}
            ),
            html.Hr(className='my-3'),
            
            # Ethnicity Distribution
            html.Label("민족 구성 (%)", className='fw-semibold text-muted mb-2'),
            dbc.Row([
                dbc.Col([
                    html.Small("East Asian", className='text-muted'),
                    dbc.Input(id='eth-asian', type='number', value=50, min=0, max=100, 
                             size='sm', className='text-center')
                ], width=4),
                dbc.Col([
                    html.Small("European", className='text-muted'),
                    dbc.Input(id='eth-european', type='number', value=30, min=0, max=100,
                             size='sm', className='text-center')
                ], width=4),
                dbc.Col([
                    html.Small("African", className='text-muted'),
                    dbc.Input(id='eth-african', type='number', value=20, min=0, max=100,
                             size='sm', className='text-center')
                ], width=4),
            ], className='mb-3'),
            
            html.Hr(className='my-3'),
            
            # Age Range
            html.Label("나이 범위", className='fw-semibold text-muted mb-1'),
            dcc.RangeSlider(
                id='age-range',
                min=18, max=80, step=1, value=[18, 65],
                marks={18: '18', 30: '30', 45: '45', 60: '60', 80: '80'},
                tooltip={'placement': 'bottom', 'always_visible': True}
            ),
            html.Hr(className='my-3'),
            
            # Gender Ratio
            html.Label("성비 (남성 %)", className='fw-semibold text-muted mb-1'),
            dcc.Slider(
                id='gender-ratio',
                min=0, max=100, step=5, value=50,
                marks={0: '0%', 25: '25%', 50: '50%', 75: '75%', 100: '100%'},
                tooltip={'placement': 'bottom', 'always_visible': True}
            ),
        ])
    ], className='shadow-sm h-100', style={'borderRadius': '12px', 'border': 'none'})


def create_physiology_card():
    """생리학적 파라미터 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("⚗️ 생리학적 파라미터", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    html.Label("체중 평균 (kg)", className='fw-semibold text-muted'),
                    dbc.Input(id='weight-mean', type='number', value=70, min=40, max=120,
                             className='mb-2')
                ], width=6),
                dbc.Col([
                    html.Label("체중 SD (kg)", className='fw-semibold text-muted'),
                    dbc.Input(id='weight-sd', type='number', value=15, min=5, max=30,
                             className='mb-2')
                ], width=6),
            ]),
            html.Hr(className='my-3'),
            html.Label("기저 간 청소율 CLint (L/h)", className='fw-semibold text-muted'),
            dbc.Input(id='base-clint', type='number', value=10, min=1, max=100,
                     className='mb-2'),
            html.Small("※ CYP2C19/3A4 유전자형에 따라 개인별로 조정됩니다", 
                      className='text-muted fst-italic')
        ])
    ], className='shadow-sm h-100', style={'borderRadius': '12px', 'border': 'none'})


def create_drug_card():
    """약물 파라미터 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("💊 약물 파라미터", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            # Drug Name with PubChem fetch
            html.Label("약물명", className='fw-semibold text-muted'),
            dbc.InputGroup([
                dbc.Input(id='drug-name', type='text', placeholder='예: Omeprazole'),
                dbc.Button("🔍 PubChem", id='fetch-pubchem', color='info', size='sm')
            ], className='mb-2'),
            html.Div(id='pubchem-status', className='mb-2'),
            
            html.Hr(className='my-3'),
            
            dbc.Row([
                dbc.Col([
                    html.Label("LogP", className='fw-semibold text-muted'),
                    dbc.Input(id='log-p', type='number', value=2.0, step=0.1)
                ], width=6),
                dbc.Col([
                    html.Label("f_u (비결합률)", className='fw-semibold text-muted'),
                    dbc.Input(id='f-u', type='number', value=0.1, min=0, max=1, step=0.01)
                ], width=6),
            ], className='mb-2'),
            
            dbc.Row([
                dbc.Col([
                    html.Label("Vd (L/kg)", className='fw-semibold text-muted'),
                    dbc.Input(id='v-d', type='number', value=1.0, step=0.1)
                ], width=6),
                dbc.Col([
                    html.Label("ka (1/h)", className='fw-semibold text-muted'),
                    dbc.Input(id='k-a', type='number', value=1.0, step=0.1)
                ], width=6),
            ], className='mb-2'),
            
            dbc.Row([
                dbc.Col([
                    html.Label("투여량 (mg)", className='fw-semibold text-muted'),
                    dbc.Input(id='dose', type='number', value=100, min=1)
                ], width=6),
                dbc.Col([
                    html.Label("생체이용률 (F)", className='fw-semibold text-muted'),
                    dbc.Input(id='bioavail', type='number', value=0.8, min=0, max=1, step=0.05)
                ], width=6),
            ]),
        ])
    ], className='shadow-sm h-100', style={'borderRadius': '12px', 'border': 'none'})


def create_simulation_control():
    """2단계 워크플로우 제어 버튼"""
    return dbc.Card([
        dbc.CardBody([
            # Step 1: 인구집단 생성 버튼
            html.Label("Step 1: 가상 인구집단 생성", className='fw-semibold text-muted mb-2'),
            dbc.Button(
                [html.Span("👥 ", style={'marginRight': '8px'}), "인구집단 생성"],
                id='generate-population',
                color='primary',
                size='lg',
                className='w-100 fw-bold mb-2',
                style={
                    'background': f'linear-gradient(135deg, {COLORS["primary"]} 0%, #2C5282 100%)',
                    'border': 'none',
                    'borderRadius': '8px',
                    'padding': '12px 24px',
                    'boxShadow': '0 4px 15px rgba(30, 58, 95, 0.4)'
                }
            ),
            dcc.Loading(
                id='loading-population',
                type='circle',
                children=html.Div(id='population-status', className='mb-3 text-center')
            ),
            
            html.Hr(className='my-3'),
            
            # Step 2: 시뮬레이션 실행 버튼
            html.Label("Step 2: 약물 시뮬레이션 실행", className='fw-semibold text-muted mb-2'),
            dbc.Button(
                [html.Span("💊 ", style={'marginRight': '8px'}), "시뮬레이션 실행"],
                id='run-simulation',
                color='success',
                size='lg',
                className='w-100 fw-bold',
                disabled=True,  # 인구집단 생성 전까지 비활성화
                style={
                    'background': f'linear-gradient(135deg, {COLORS["secondary"]} 0%, #27AE60 100%)',
                    'border': 'none',
                    'borderRadius': '8px',
                    'padding': '12px 24px',
                    'boxShadow': '0 4px 15px rgba(46, 204, 113, 0.4)'
                }
            ),
            dcc.Loading(
                id='loading-simulation',
                type='circle',
                children=html.Div(id='simulation-status', className='mt-2 text-center')
            )
        ])
    ], className='shadow-sm', style={'borderRadius': '12px', 'border': 'none'})


def create_pk_curves_card():
    """PK 곡선 시각화 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("📈 혈장 농도-시간 곡선", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            dcc.Graph(
                id='pk-curves',
                config={'displayModeBar': True, 'responsive': True},
                style={'height': '400px'}
            )
        ])
    ], className='shadow-sm', style={'borderRadius': '12px', 'border': 'none'})


def create_safety_card():
    """안전 마진 분석 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("⚠️ 안전 마진 분석", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            html.Label("독성 임계값 (ng/mL)", className='fw-semibold text-muted'),
            dbc.Input(id='toxic-threshold', type='number', value=1000, min=1,
                     className='mb-3'),
            html.Div(id='safety-report', className='mt-3'),
            dcc.Graph(
                id='cmax-histogram',
                config={'displayModeBar': False, 'responsive': True},
                style={'height': '250px'}
            )
        ])
    ], className='shadow-sm', style={'borderRadius': '12px', 'border': 'none'})


def create_population_summary_card():
    """인구집단 요약 카드"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("📊 인구집단 요약", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'})
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            html.Div(id='population-summary'),
            dcc.Graph(
                id='metabolizer-pie',
                config={'displayModeBar': False, 'responsive': True},
                style={'height': '250px'}
            )
        ])
    ], className='shadow-sm', style={'borderRadius': '12px', 'border': 'none'})


def create_digital_twin_card():
    """디지털 트윈 시각화 카드 - 사람 모양 아이콘으로 인구집단 표시"""
    return dbc.Card([
        dbc.CardHeader([
            html.H5("🧬 디지털 트윈 인구집단", className='mb-0',
                   style={'color': COLORS['primary'], 'fontWeight': '600'}),
        ], style={'backgroundColor': '#F1F5F9'}),
        dbc.CardBody([
            # 범례
            html.Div([
                html.Span("대사자 표현형: ", className='fw-semibold me-2'),
                html.Span("🟢 NM", className='me-2', style={'fontSize': '0.85rem'}),
                html.Span("🟡 IM", className='me-2', style={'fontSize': '0.85rem'}),
                html.Span("🔴 PM", className='me-2', style={'fontSize': '0.85rem'}),
                html.Span("🔵 UM", style={'fontSize': '0.85rem'}),
            ], className='mb-3 text-muted', style={'fontSize': '0.9rem'}),
            
            # 디지털 트윈 그리드 컨테이너
            html.Div(
                id='digital-twin-grid',
                style={
                    'display': 'flex',
                    'flexWrap': 'wrap',
                    'gap': '4px',
                    'justifyContent': 'center',
                    'maxHeight': '350px',
                    'overflowY': 'auto',
                    'padding': '10px',
                    'backgroundColor': '#F8FAFC',
                    'borderRadius': '8px'
                }
            ),
            
            # 선택된 개인 정보 표시
            html.Div(id='selected-individual-info', className='mt-3')
        ])
    ], className='shadow-sm', style={'borderRadius': '12px', 'border': 'none'})


def create_person_icon(individual_data: dict, index: int) -> html.Div:
    """개인을 나타내는 사람 모양 아이콘 생성
    
    Args:
        individual_data: 개인 정보 딕셔너리
        index: 개인 인덱스
        
    Returns:
        html.Div 컴포넌트
    """
    # 대사자 표현형에 따른 색상
    metabolizer = individual_data.get('metabolizer', 'NM')
    gender = individual_data.get('gender', 'M')
    
    color_map = {
        'PM': '#E74C3C',   # 빨강
        'IM': '#F39C12',   # 노랑
        'NM': '#2ECC71',   # 초록
        'UM': '#3498DB'    # 파랑
    }
    color = color_map.get(metabolizer, '#2ECC71')
    
    # 성별에 따른 아이콘 (SVG)
    if gender == 'M':
        icon = '👨'
    else:
        icon = '👩'
    
    # 툴팁 텍스트
    age = individual_data.get('age', 0)
    weight = individual_data.get('weight', 0)
    ethnicity = individual_data.get('ethnicity', '')
    
    tooltip_text = f"ID: {index+1}\n나이: {age}세\n체중: {weight:.1f}kg\n민족: {ethnicity}\n표현형: {metabolizer}"
    
    return html.Div(
        icon,
        id={'type': 'person-icon', 'index': index},
        title=tooltip_text,
        style={
            'width': '28px',
            'height': '28px',
            'display': 'flex',
            'alignItems': 'center',
            'justifyContent': 'center',
            'fontSize': '18px',
            'backgroundColor': color,
            'borderRadius': '50%',
            'cursor': 'pointer',
            'transition': 'transform 0.2s, box-shadow 0.2s',
            'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
        },
        className='person-icon'
    )


# ============================================================================
# Main Layout
# ============================================================================

app.layout = html.Div([
    create_header(),
    
    dbc.Container([
        dbc.Row([
            # Left Panel - Input Controls
            dbc.Col([
                create_population_card(),
                html.Div(className='mb-3'),
                create_physiology_card(),
                html.Div(className='mb-3'),
                create_drug_card(),
                html.Div(className='mb-3'),
                create_simulation_control(),
            ], lg=4, md=5, className='mb-4'),
            
            # Right Panel - Visualizations
            dbc.Col([
                # 디지털 트윈 시각화 (상단)
                create_digital_twin_card(),
                html.Div(className='mb-3'),
                
                # PK 곡선 (중단)
                create_pk_curves_card(),
                html.Div(className='mb-3'),
                
                # 하단 패널들
                dbc.Row([
                    dbc.Col([create_safety_card()], lg=6),
                    dbc.Col([create_population_summary_card()], lg=6),
                ])
            ], lg=8, md=7),
        ])
    ], fluid=True, style={'maxWidth': '1600px'}),
    
    # Store for simulation results
    dcc.Store(id='simulation-results'),
    dcc.Store(id='population-individuals'),  # 개인별 데이터 저장
    
], style={'backgroundColor': COLORS['background'], 'minHeight': '100vh', 'paddingBottom': '50px'})


# ============================================================================
# Callbacks
# ============================================================================

@callback(
    [Output('pubchem-status', 'children'),
     Output('log-p', 'value')],
    Input('fetch-pubchem', 'n_clicks'),
    State('drug-name', 'value'),
    prevent_initial_call=True
)
def fetch_drug_data(n_clicks, drug_name):
    """PubChem에서 약물 데이터 가져오기"""
    if not drug_name:
        return dbc.Alert("약물명을 입력하세요", color='warning', className='py-1 mb-0'), 2.0
    
    data = fetch_pubchem_data(drug_name)
    
    if data['found']:
        return (
            dbc.Alert(
                f"✓ {drug_name} 발견 (MW: {data['mw']:.1f})",
                color='success', className='py-1 mb-0'
            ),
            round(data['log_p'], 2)
        )
    else:
        return dbc.Alert("약물을 찾을 수 없습니다", color='danger', className='py-1 mb-0'), 2.0


@callback(
    [Output('population-individuals', 'data'),
     Output('population-status', 'children'),
     Output('run-simulation', 'disabled')],
    Input('generate-population', 'n_clicks'),
    [State('n-subjects', 'value'),
     State('eth-asian', 'value'),
     State('eth-european', 'value'),
     State('eth-african', 'value'),
     State('age-range', 'value'),
     State('gender-ratio', 'value'),
     State('weight-mean', 'value'),
     State('weight-sd', 'value'),
     State('base-clint', 'value')],
    prevent_initial_call=True
)
def generate_population(n_clicks, n_subjects, eth_asian, eth_european, eth_african,
                        age_range, gender_ratio, weight_mean, weight_sd, base_clint):
    """Step 1: 가상 인구집단 생성"""
    
    # Normalize ethnicity percentages
    total_eth = eth_asian + eth_european + eth_african
    if total_eth == 0:
        total_eth = 100
    
    ethnicity_dist = {
        Ethnicity.EAST_ASIAN: eth_asian / total_eth,
        Ethnicity.EUROPEAN: eth_european / total_eth,
        Ethnicity.AFRICAN: eth_african / total_eth
    }
    
    # Generate population
    generator = PopulationGenerator(
        n_subjects=n_subjects,
        ethnicity_distribution=ethnicity_dist,
        age_range=tuple(age_range),
        gender_ratio=gender_ratio / 100,
        weight_mean=weight_mean,
        weight_sd=weight_sd,
        base_cl_int=base_clint,
        random_seed=None  # 매번 다른 인구집단 생성
    )
    
    population = generator.generate()
    pop_summary = generator.get_population_summary(population)
    
    # 개인별 데이터 저장 (디지털 트윈 + 시뮬레이션용)
    individuals_data = []
    for ind in population:
        # 대사자 표현형 약여
        metabolizer_map = {
            'Poor Metabolizer (PM)': 'PM',
            'Intermediate Metabolizer (IM)': 'IM',
            'Normal Metabolizer (NM)': 'NM',
            'Ultra-rapid Metabolizer (UM)': 'UM'
        }
        individuals_data.append({
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
            # 생리학적 파라미터도 저장 (시뮬레이션용)
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
    
    # pop_summary도 함께 저장
    result_data = {
        'individuals': individuals_data,
        'summary': pop_summary
    }
    
    status = dbc.Alert(
        f"✓ {n_subjects}명 가상 인구집단 생성 완료! 약물 정보를 입력하고 시뮬레이션을 실행하세요.",
        color='success', className='py-1 mb-0'
    )
    
    # 시뮬레이션 버튼 활성화 (disabled=False)
    return result_data, status, False


@callback(
    [Output('simulation-results', 'data'),
     Output('simulation-status', 'children')],
    Input('run-simulation', 'n_clicks'),
    [State('population-individuals', 'data'),
     State('drug-name', 'value'),
     State('log-p', 'value'),
     State('f-u', 'value'),
     State('v-d', 'value'),
     State('k-a', 'value'),
     State('dose', 'value'),
     State('bioavail', 'value')],
    prevent_initial_call=True
)
def run_drug_simulation(n_clicks, population_data, drug_name, log_p, f_u, v_d, k_a, dose, bioavail):
    """Step 2: 약물 시뮬레이션 실행"""
    
    if population_data is None:
        return None, dbc.Alert("먼저 인구집단을 생성하세요!", color='warning', className='py-1 mb-0')
    
    individuals = population_data.get('individuals', [])
    pop_summary = population_data.get('summary', {})
    
    if not individuals:
        return None, dbc.Alert("인구집단 데이터가 없습니다.", color='danger', className='py-1 mb-0')
    
    # Drug parameters
    drug_params = DrugParameters(
        name=drug_name or "Unknown Drug",
        log_p=log_p,
        f_u=f_u,
        v_d=v_d,
        k_a=k_a,
        f=bioavail
    )
    
    # Simulation config
    sim_config = SimulationConfig(
        dose=dose,
        route="oral",
        t_max=24,
        n_points=241
    )
    
    # 저장된 생리학적 파라미터로 PhysiologicalParameters 객체 생성
    population_phys = []
    for ind in individuals:
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
    
    # Run population simulation
    sim_results = run_population_simulation(drug_params, population_phys, sim_config)
    
    # Prepare data for storage
    results = {
        'time': sim_results['time'].tolist(),
        'mean_concentration': sim_results['mean_concentration'].tolist(),
        'ci_lower': sim_results['ci_lower'].tolist(),
        'ci_upper': sim_results['ci_upper'].tolist(),
        'individual_curves': sim_results['individual_curves'][:50].tolist(),
        'cmax_distribution': sim_results['cmax_distribution'].tolist(),
        'auc_distribution': sim_results['auc_distribution'].tolist(),
        'pop_summary': pop_summary
    }
    
    status = dbc.Alert(
        f"✓ {len(individuals)}명 대상 '{drug_name or 'Drug'}' 시뮬레이션 완료!",
        color='success', className='py-1 mb-0'
    )
    
    return results, status


@callback(
    Output('pk-curves', 'figure'),
    Input('simulation-results', 'data')
)
def update_pk_curves(results):
    """PK 곡선 업데이트"""
    
    if results is None:
        # Empty plot
        fig = go.Figure()
        fig.add_annotation(
            text="시뮬레이션을 실행하세요",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16, color=COLORS['text_muted'])
        )
        fig.update_layout(
            template='plotly_white',
            margin=dict(l=40, r=40, t=40, b=40)
        )
        return fig
    
    time = results['time']
    mean_conc = results['mean_concentration']
    ci_lower = results['ci_lower']
    ci_upper = results['ci_upper']
    individual_curves = results['individual_curves']
    
    fig = go.Figure()
    
    # Individual curves (subset)
    for i, curve in enumerate(individual_curves[:30]):
        fig.add_trace(go.Scatter(
            x=time, y=curve,
            mode='lines',
            line=dict(color='rgba(150, 150, 150, 0.2)', width=0.5),
            name='Individual' if i == 0 else None,
            showlegend=(i == 0),
            hoverinfo='skip'
        ))
    
    # 90% CI band
    fig.add_trace(go.Scatter(
        x=time + time[::-1],
        y=ci_upper + ci_lower[::-1],
        fill='toself',
        fillcolor='rgba(30, 58, 95, 0.2)',
        line=dict(color='rgba(0,0,0,0)'),
        name='90% CI',
        hoverinfo='skip'
    ))
    
    # Mean curve
    fig.add_trace(go.Scatter(
        x=time, y=mean_conc,
        mode='lines',
        line=dict(color=COLORS['primary'], width=3),
        name='집단 평균'
    ))
    
    fig.update_layout(
        template='plotly_white',
        xaxis_title='시간 (h)',
        yaxis_title='혈장 농도 (ng/mL)',
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='right',
            x=1
        ),
        margin=dict(l=60, r=40, t=60, b=60),
        hovermode='x unified'
    )
    
    return fig


@callback(
    [Output('safety-report', 'children'),
     Output('cmax-histogram', 'figure')],
    [Input('simulation-results', 'data'),
     Input('toxic-threshold', 'value')]
)
def update_safety_analysis(results, toxic_threshold):
    """안전 마진 분석 업데이트"""
    
    empty_fig = go.Figure()
    empty_fig.update_layout(template='plotly_white', margin=dict(l=40, r=40, t=20, b=40))
    
    if results is None:
        return html.Div("시뮬레이션 결과 대기 중...", className='text-muted'), empty_fig
    
    cmax_dist = np.array(results['cmax_distribution'])
    safety = calculate_safety_margin(cmax_dist, toxic_threshold)
    
    # Safety report
    if safety['percentage_exceeding'] > 10:
        alert_color = 'danger'
        icon = '🔴'
    elif safety['percentage_exceeding'] > 5:
        alert_color = 'warning'
        icon = '🟡'
    else:
        alert_color = 'success'
        icon = '🟢'
    
    report = dbc.Alert([
        html.H6(f"{icon} 안전 마진 분석 결과", className='alert-heading'),
        html.Hr(className='my-2'),
        html.P([
            html.Strong(f"{safety['percentage_exceeding']:.1f}%"),
            f" 인구가 독성 임계값 ({toxic_threshold} ng/mL) 초과"
        ], className='mb-1'),
        html.P([
            "95th percentile Cmax: ",
            html.Strong(f"{safety['cmax_95th_percentile']:.1f} ng/mL")
        ], className='mb-1'),
        html.P([
            "안전 비율 (Threshold/95th): ",
            html.Strong(f"{safety['safety_ratio']:.2f}x")
        ], className='mb-0')
    ], color=alert_color, className='mb-0')
    
    # Histogram
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=cmax_dist,
        nbinsx=30,
        marker_color=COLORS['primary'],
        opacity=0.7,
        name='Cmax 분포'
    ))
    
    # Toxic threshold line
    fig.add_vline(
        x=toxic_threshold,
        line_dash='dash',
        line_color=COLORS['danger'],
        annotation_text='독성 임계값',
        annotation_position='top right'
    )
    
    fig.update_layout(
        template='plotly_white',
        xaxis_title='Cmax (ng/mL)',
        yaxis_title='빈도',
        margin=dict(l=40, r=20, t=30, b=40),
        showlegend=False
    )
    
    return report, fig


@callback(
    [Output('population-summary', 'children'),
     Output('metabolizer-pie', 'figure')],
    Input('simulation-results', 'data')
)
def update_population_summary(results):
    """인구집단 요약 업데이트"""
    
    empty_fig = go.Figure()
    empty_fig.update_layout(template='plotly_white', margin=dict(l=20, r=20, t=20, b=20))
    
    if results is None:
        return html.Div("시뮬레이션 결과 대기 중...", className='text-muted'), empty_fig
    
    summary = results['pop_summary']
    
    # Summary badges
    summary_div = html.Div([
        dbc.Badge(f"N = {summary['n_subjects']}", color='primary', className='me-1 mb-1'),
        dbc.Badge(
            f"나이: {summary['demographics']['age']['mean']:.0f} ± {summary['demographics']['age']['sd']:.0f}",
            color='info', className='me-1 mb-1'
        ),
        dbc.Badge(
            f"체중: {summary['demographics']['weight']['mean']:.0f} ± {summary['demographics']['weight']['sd']:.0f} kg",
            color='info', className='me-1 mb-1'
        ),
        dbc.Badge(
            f"남성: {summary['demographics']['gender']['male_ratio']*100:.0f}%",
            color='secondary', className='mb-1'
        ),
    ])
    
    # Metabolizer pie chart
    met_dist = summary['metabolizer_distribution']
    labels = ['PM', 'IM', 'NM', 'UM']
    values = [
        met_dist.get('Poor Metabolizer (PM)', 0),
        met_dist.get('Intermediate Metabolizer (IM)', 0),
        met_dist.get('Normal Metabolizer (NM)', 0),
        met_dist.get('Ultra-rapid Metabolizer (UM)', 0)
    ]
    
    colors = [COLORS['danger'], COLORS['warning'], COLORS['secondary'], COLORS['primary']]
    
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.4,
        marker_colors=colors,
        textinfo='label+percent',
        textposition='outside'
    )])
    
    fig.update_layout(
        template='plotly_white',
        margin=dict(l=20, r=20, t=30, b=20),
        showlegend=False,
        annotations=[dict(text='대사자<br>표현형', x=0.5, y=0.5, font_size=11, showarrow=False)]
    )
    
    return summary_div, fig


@callback(
    Output('digital-twin-grid', 'children'),
    Input('population-individuals', 'data')
)
def update_digital_twin_grid(population_data):
    """디지털 트윈 그리드 업데이트 - 사람 모양 아이콘으로 인구집단 시각화"""
    
    if population_data is None:
        return html.Div([
            html.P("👆 Step 1: '인구집단 생성' 버튼을 클릭하세요",
                  className='text-muted text-center py-5',
                  style={'fontSize': '1.1rem'})
        ])
    
    # 새 데이터 구조에서 individuals 추출
    individuals = population_data.get('individuals', [])
    
    if not individuals:
        return html.Div([
            html.P("인구집단 데이터가 없습니다.", className='text-muted text-center py-5')
        ])
    
    # 최대 200명까지만 표시 (성능 고려)
    display_limit = min(200, len(individuals))
    
    icons = []
    for i, ind in enumerate(individuals[:display_limit]):
        # 대사자 표현형에 따른 색상
        metabolizer = ind.get('metabolizer', 'NM')
        gender = ind.get('gender', 'M')
        
        color_map = {
            'PM': '#E74C3C',   # 빨강
            'IM': '#F39C12',   # 노랑
            'NM': '#2ECC71',   # 초록
            'UM': '#3498DB'    # 파랑
        }
        color = color_map.get(metabolizer, '#2ECC71')
        
        # 성별에 따른 아이콘
        icon = '👨' if gender == 'M' else '👩'
        
        # 툴팁
        tooltip = f"ID: {ind['id']} | {ind['age']}세 | {ind['weight']}kg | {ind['ethnicity']} | {metabolizer}"
        
        icons.append(
            html.Div(
                icon,
                title=tooltip,
                style={
                    'width': '26px',
                    'height': '26px',
                    'display': 'flex',
                    'alignItems': 'center',
                    'justifyContent': 'center',
                    'fontSize': '14px',
                    'backgroundColor': color,
                    'borderRadius': '50%',
                    'cursor': 'pointer',
                    'transition': 'transform 0.2s',
                    'boxShadow': '0 1px 3px rgba(0,0,0,0.15)'
                }
            )
        )
    
    # 더 많은 개인이 있을 경우 표시
    if len(individuals) > display_limit:
        icons.append(
            html.Div(
                f"+{len(individuals) - display_limit}",
                style={
                    'padding': '4px 10px',
                    'fontSize': '12px',
                    'color': COLORS['text_muted'],
                    'fontWeight': '600'
                }
            )
        )
    
    return icons


# ============================================================================
# Run Server
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*60)
    print("  🧬 Virtual Population PBPK Simulator")
    print("="*60)
    print("  서버 시작 중...")
    print("  브라우저에서 http://localhost:8050 을 열어주세요")
    print("="*60 + "\n")
    
    app.run(debug=True, port=8050)
