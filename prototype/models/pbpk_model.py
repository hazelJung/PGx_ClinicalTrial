"""
PBPK Model - Multi-compartment Physiologically Based Pharmacokinetic Model
============================================================================
3-구획 PBPK 모델 (장관, 중심, 간)을 구현합니다.

Compartments:
    1. Gut (장관) - 경구 투여 시 약물 흡수
    2. Central/Plasma (중심/혈장) - 전신 순환
    3. Liver (간) - 대사적 청소율
"""

import numpy as np
from scipy.integrate import odeint, trapezoid
from dataclasses import dataclass, field
from typing import Tuple, Optional, Dict, Any


@dataclass
class DrugParameters:
    """약물 파라미터 데이터 클래스
    
    Attributes:
        name: 약물명
        log_p: 지질/물 분배계수 (LogP)
        f_u: 혈장 단백 비결합률 (fraction unbound)
        v_d: 분포용적 (L/kg)
        k_a: 흡수속도상수 (1/h)
        f: 생체이용률 (bioavailability)
        mw: 분자량 (g/mol)
    """
    name: str = "Generic Drug"
    log_p: float = 2.0
    f_u: float = 0.1
    v_d: float = 1.0  # L/kg
    k_a: float = 1.0  # 1/h
    f: float = 0.8    # bioavailability
    mw: float = 300.0 # molecular weight


@dataclass
class PhysiologicalParameters:
    """생리학적 파라미터 데이터 클래스
    
    Attributes:
        body_weight: 체중 (kg)
        v_plasma: 혈장 용적 (L)
        v_liver: 간 용적 (L)
        q_liver: 간 혈류량 (L/h)
        cl_int: 간 내재적 청소율 (L/h) - 기저 값
        cl_renal: 신장 청소율 (L/h)
        activity_score: CYP 효소 활성 점수 (0, 0.5, 1.0, 2.0)
    """
    body_weight: float = 70.0  # kg
    v_plasma: float = 3.0      # L
    v_liver: float = 1.5       # L
    q_liver: float = 90.0      # L/h (hepatic blood flow)
    cl_int: float = 10.0       # L/h (intrinsic clearance)
    cl_renal: float = 0.0      # L/h
    activity_score: float = 1.0  # metabolizer status


@dataclass
class SimulationConfig:
    """시뮬레이션 설정
    
    Attributes:
        dose: 투여 용량 (mg)
        route: 투여 경로 ('oral' or 'iv')
        t_max: 시뮬레이션 종료 시간 (h)
        n_points: 시간 포인트 수
    """
    dose: float = 100.0       # mg
    route: str = "oral"       # 'oral' or 'iv'
    t_max: float = 24.0       # hours
    n_points: int = 241       # time points


class PBPKModel:
    """3-구획 PBPK 모델
    
    Multi-compartment PBPK model for drug PK simulation.
    Implements Gut → Plasma → Liver compartmental flow.
    """
    
    def __init__(
        self,
        drug_params: Optional[DrugParameters] = None,
        phys_params: Optional[PhysiologicalParameters] = None,
        sim_config: Optional[SimulationConfig] = None
    ):
        self.drug = drug_params or DrugParameters()
        self.phys = phys_params or PhysiologicalParameters()
        self.config = sim_config or SimulationConfig()
        
        # Liver-plasma partition coefficient (Kp) estimation using LogP
        # Simplified Poulin-Theil method
        self.k_p_liver = self._estimate_kp_liver()
        
    def _estimate_kp_liver(self) -> float:
        """LogP를 기반으로 간-혈장 분배계수(Kp) 추정
        
        Simplified estimation based on lipophilicity.
        """
        # Poulin-Theil simplified approach
        log_p = self.drug.log_p
        f_u = self.drug.f_u
        
        # Neutral lipophilic drug partition
        kp = 0.5 + 0.5 * 10 ** (0.7 * log_p - 0.3) * f_u
        return max(1.0, min(kp, 50.0))  # Clamp between 1-50
    
    def _ode_system(self, y: np.ndarray, t: float, params: Dict[str, float]) -> list:
        """PBPK ODE 시스템
        
        State variables:
            y[0]: A_gut - 장관 내 약물량 (mg)
            y[1]: C_plasma - 혈장 농도 (mg/L = μg/mL)
            y[2]: C_liver - 간 농도 (mg/L)
        """
        A_gut, C_plasma, C_liver = y
        
        # Parameters
        k_a = params['k_a']
        F = params['F']
        V_c = params['V_c']
        V_liver = params['V_liver']
        Q_liver = params['Q_liver']
        CL_int = params['CL_int']
        CL_renal = params['CL_renal']
        f_u = params['f_u']
        K_p = params['K_p']
        
        # Differential equations
        # 장관 구획: 1차 흡수
        dA_gut_dt = -k_a * A_gut
        
        # 중심 구획 (혈장)
        absorption = (k_a * A_gut * F) / V_c
        liver_uptake = (Q_liver / V_c) * C_plasma
        liver_return = (Q_liver / V_c) * (C_liver / K_p)
        renal_clearance = (CL_renal / V_c) * C_plasma
        
        dC_plasma_dt = absorption - liver_uptake + liver_return - renal_clearance
        
        # 간 구획
        hepatic_uptake = (Q_liver / V_liver) * C_plasma
        hepatic_return = (Q_liver / V_liver) * (C_liver / K_p)
        metabolism = (CL_int * f_u / V_liver) * C_liver
        
        dC_liver_dt = hepatic_uptake - hepatic_return - metabolism
        
        return [dA_gut_dt, dC_plasma_dt, dC_liver_dt]
    
    def solve(self) -> Dict[str, Any]:
        """PBPK 모델 풀이
        
        Returns:
            Dictionary containing:
                - time: 시간 배열 (h)
                - c_plasma: 혈장 농도 배열 (ng/mL)
                - c_liver: 간 농도 배열 (ng/mL)
                - pk_metrics: PK 파라미터 (Cmax, Tmax, AUC, t_half)
        """
        # Time array
        t = np.linspace(0, self.config.t_max, self.config.n_points)
        
        # Effective intrinsic clearance (adjusted by activity score)
        effective_cl_int = self.phys.cl_int * self.phys.activity_score
        
        # Central volume (adjusted for body weight)
        v_central = self.drug.v_d * self.phys.body_weight
        
        # Parameter dictionary for ODE
        params = {
            'k_a': self.drug.k_a,
            'F': self.drug.f,
            'V_c': v_central,
            'V_liver': self.phys.v_liver,
            'Q_liver': self.phys.q_liver,
            'CL_int': effective_cl_int,
            'CL_renal': self.phys.cl_renal,
            'f_u': self.drug.f_u,
            'K_p': self.k_p_liver
        }
        
        # Initial conditions
        if self.config.route == "oral":
            # Oral: all drug in gut
            y0 = [self.config.dose, 0.0, 0.0]
        else:
            # IV: all drug in plasma
            y0 = [0.0, self.config.dose / v_central, 0.0]
        
        # Solve ODE
        solution = odeint(self._ode_system, y0, t, args=(params,))
        
        # Extract concentrations (convert to ng/mL: mg/L * 1000 = μg/L = ng/mL)
        c_plasma = solution[:, 1] * 1000  # ng/mL
        c_liver = solution[:, 2] * 1000   # ng/mL
        
        # Calculate PK metrics
        pk_metrics = self._calculate_pk_metrics(t, c_plasma)
        
        return {
            'time': t,
            'c_plasma': c_plasma,
            'c_liver': c_liver,
            'pk_metrics': pk_metrics,
            'parameters': params
        }
    
    def _calculate_pk_metrics(self, t: np.ndarray, c_plasma: np.ndarray) -> Dict[str, float]:
        """PK 파라미터 계산
        
        Args:
            t: 시간 배열
            c_plasma: 혈장 농도 배열 (ng/mL)
            
        Returns:
            Dictionary with Cmax, Tmax, AUC, t_half
        """
        # Cmax and Tmax
        cmax_idx = np.argmax(c_plasma)
        cmax = float(c_plasma[cmax_idx])
        tmax = float(t[cmax_idx])
        
        # AUC (trapezoidal rule)
        auc = float(trapezoid(c_plasma, t))
        
        # Half-life estimation (terminal phase)
        # Find points after Cmax where concentration > 10% of Cmax
        terminal_mask = (t > tmax) & (c_plasma > 0.1 * cmax)
        if np.sum(terminal_mask) > 2:
            t_terminal = t[terminal_mask]
            c_terminal = c_plasma[terminal_mask]
            # Log-linear regression
            log_c = np.log(c_terminal + 1e-10)
            slope, _ = np.polyfit(t_terminal, log_c, 1)
            t_half = -np.log(2) / slope if slope < 0 else np.nan
        else:
            t_half = np.nan
        
        return {
            'cmax': cmax,          # ng/mL
            'tmax': tmax,          # h
            'auc': auc,            # ng·h/mL
            't_half': float(t_half) if not np.isnan(t_half) else None  # h
        }


def run_population_simulation(
    drug_params: DrugParameters,
    population_phys_params: list,
    sim_config: SimulationConfig
) -> Dict[str, Any]:
    """집단 시뮬레이션 실행
    
    Args:
        drug_params: 약물 파라미터
        population_phys_params: 개인별 생리학적 파라미터 리스트
        sim_config: 시뮬레이션 설정
        
    Returns:
        Dictionary with population simulation results
    """
    n_subjects = len(population_phys_params)
    n_points = sim_config.n_points
    
    # Result arrays
    all_c_plasma = np.zeros((n_subjects, n_points))
    all_pk_metrics = []
    
    time = None
    
    for i, phys_params in enumerate(population_phys_params):
        model = PBPKModel(drug_params, phys_params, sim_config)
        result = model.solve()
        
        if time is None:
            time = result['time']
        
        all_c_plasma[i, :] = result['c_plasma']
        all_pk_metrics.append(result['pk_metrics'])
    
    # Population statistics
    mean_c = np.mean(all_c_plasma, axis=0)
    std_c = np.std(all_c_plasma, axis=0)
    percentile_5 = np.percentile(all_c_plasma, 5, axis=0)
    percentile_95 = np.percentile(all_c_plasma, 95, axis=0)
    
    # Extract Cmax values for safety analysis
    cmax_values = [m['cmax'] for m in all_pk_metrics]
    auc_values = [m['auc'] for m in all_pk_metrics]
    
    return {
        'time': time,
        'individual_curves': all_c_plasma,
        'mean_concentration': mean_c,
        'std_concentration': std_c,
        'ci_lower': percentile_5,
        'ci_upper': percentile_95,
        'cmax_distribution': np.array(cmax_values),
        'auc_distribution': np.array(auc_values),
        'pk_metrics_list': all_pk_metrics
    }


if __name__ == "__main__":
    # Test the model
    print("=== PBPK Model Test ===\n")
    
    drug = DrugParameters(
        name="Test Drug",
        log_p=2.5,
        f_u=0.15,
        v_d=1.2,
        k_a=0.8,
        f=0.75
    )
    
    phys = PhysiologicalParameters(
        body_weight=70,
        activity_score=1.0
    )
    
    config = SimulationConfig(
        dose=100,
        route="oral",
        t_max=24,
        n_points=241
    )
    
    model = PBPKModel(drug, phys, config)
    result = model.solve()
    
    print(f"Drug: {drug.name}")
    print(f"Dose: {config.dose} mg ({config.route})")
    print(f"\nPK Parameters:")
    print(f"  Cmax: {result['pk_metrics']['cmax']:.2f} ng/mL")
    print(f"  Tmax: {result['pk_metrics']['tmax']:.2f} h")
    print(f"  AUC: {result['pk_metrics']['auc']:.2f} ng·h/mL")
    if result['pk_metrics']['t_half']:
        print(f"  t1/2: {result['pk_metrics']['t_half']:.2f} h")
