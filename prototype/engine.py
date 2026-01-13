"""
Population Engine - Virtual Population Generator with Genetic Variability
===========================================================================
가상 인구집단 생성기 (Monte Carlo 샘플링 + 약물유전체학)

Features:
    - 민족별 CYP 효소 유전자 빈도 매핑 (gnomAD/PharmGKB 기반)
    - Hardy-Weinberg 평형을 이용한 유전자형 할당
    - Activity Score 기반 대사 청소율 조정
    - 체중-장기 용적 상관관계 모델링
"""

import json
import os
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from enum import Enum
from pathlib import Path

from pbpk_model import PhysiologicalParameters


# ============================================================================
# Data Loading from JSON Files
# ============================================================================

def get_data_dir() -> Path:
    """데이터 디렉토리 경로 반환"""
    return Path(__file__).parent / "data"


def load_allele_frequencies() -> Dict:
    """JSON 파일에서 대립유전자 빈도 데이터 로드
    
    Returns:
        민족별 CYP 효소 대립유전자 빈도 딕셔너리
    """
    data_path = get_data_dir() / "allele_frequencies.json"
    
    if not data_path.exists():
        raise FileNotFoundError(
            f"Allele frequency data file not found: {data_path}\n"
            "Please run the data fetch script first."
        )
    
    with open(data_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    return data.get("allele_frequencies", {})


def load_activity_scores() -> Dict:
    """JSON 파일에서 Activity Score 데이터 로드
    
    Returns:
        CYP 효소별 대립유전자 Activity Score 딕셔너리
    """
    data_path = get_data_dir() / "activity_scores.json"
    
    if not data_path.exists():
        raise FileNotFoundError(
            f"Activity score data file not found: {data_path}\n"
            "Please run the data fetch script first."
        )
    
    with open(data_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    return data.get("activity_scores", {})


# ============================================================================
# Enums and Data Classes  
# ============================================================================

class Ethnicity(Enum):
    """민족 분류"""
    EAST_ASIAN = "East Asian"
    EUROPEAN = "European"
    AFRICAN = "African"
    LATINO = "Latino"
    CENTRAL_SOUTH_ASIAN = "Central/South Asian"


class MetabolizerStatus(Enum):
    """대사자 표현형"""
    POOR = "Poor Metabolizer (PM)"
    INTERMEDIATE = "Intermediate Metabolizer (IM)"
    NORMAL = "Normal Metabolizer (NM)"
    ULTRA_RAPID = "Ultra-rapid Metabolizer (UM)"


# Load data from JSON files at module import
# This allows the data to be cached and reused
try:
    _ALLELE_FREQ_DATA = load_allele_frequencies()
    _ACTIVITY_SCORE_DATA = load_activity_scores()
except FileNotFoundError as e:
    print(f"Warning: {e}")
    _ALLELE_FREQ_DATA = {}
    _ACTIVITY_SCORE_DATA = {}


def get_allele_frequencies(ethnicity: Ethnicity, gene: str) -> Dict[str, float]:
    """특정 민족과 유전자에 대한 대립유전자 빈도 반환
    
    Args:
        ethnicity: 민족 Enum
        gene: 유전자명 (예: 'CYP2C19')
        
    Returns:
        대립유전자별 빈도 딕셔너리
    """
    eth_name = ethnicity.value
    
    if eth_name not in _ALLELE_FREQ_DATA:
        # 데이터가 없으면 기본값으로 East Asian 사용
        eth_name = "East Asian"
    
    if gene not in _ALLELE_FREQ_DATA.get(eth_name, {}):
        # 유전자 데이터가 없으면 기본값 반환
        return {"*1": 1.0}
    
    return _ALLELE_FREQ_DATA[eth_name][gene]


def get_activity_score(gene: str, allele: str) -> float:
    """특정 유전자와 대립유전자에 대한 Activity Score 반환
    
    Args:
        gene: 유전자명 (예: 'CYP2C19')
        allele: 대립유전자명 (예: '*1')
        
    Returns:
        Activity Score (기본값: 1.0)
    """
    if gene not in _ACTIVITY_SCORE_DATA:
        return 1.0
    
    return _ACTIVITY_SCORE_DATA[gene].get(allele, 1.0)


@dataclass
class IndividualCharacteristics:
    """개인 특성 데이터 클래스"""
    subject_id: int
    age: int
    gender: str  # 'M' or 'F'
    ethnicity: Ethnicity
    weight: float  # kg
    height: float  # cm
    bmi: float
    
    # 유전체 정보
    cyp2c19_genotype: Tuple[str, str]
    cyp3a4_genotype: Tuple[str, str]
    cyp2c19_activity_score: float
    cyp3a4_activity_score: float
    combined_activity_score: float
    metabolizer_status: MetabolizerStatus
    
    # 생리학적 파라미터
    phys_params: PhysiologicalParameters = None


class PopulationGenerator:
    """가상 인구집단 생성기
    
    Monte Carlo 샘플링을 사용하여 생리학적, 유전적 다양성을 가진
    가상 인구집단을 생성합니다.
    """
    
    def __init__(
        self,
        n_subjects: int = 1000,
        ethnicity_distribution: Optional[Dict[Ethnicity, float]] = None,
        age_range: Tuple[int, int] = (18, 65),
        gender_ratio: float = 0.5,  # Male ratio
        weight_mean: float = 70.0,
        weight_sd: float = 15.0,
        base_cl_int: float = 10.0,
        random_seed: Optional[int] = None
    ):
        """
        Args:
            n_subjects: 생성할 개인 수
            ethnicity_distribution: 민족별 비율 (예: {EAST_ASIAN: 0.5, EUROPEAN: 0.3, AFRICAN: 0.2})
            age_range: 나이 범위 (min, max)
            gender_ratio: 남성 비율 (0-1)
            weight_mean: 평균 체중 (kg)
            weight_sd: 체중 표준편차 (kg)
            base_cl_int: 기저 간 내재적 청소율 (L/h)
            random_seed: 재현성을 위한 난수 시드
        """
        self.n_subjects = n_subjects
        self.age_range = age_range
        self.gender_ratio = gender_ratio
        self.weight_mean = weight_mean
        self.weight_sd = weight_sd
        self.base_cl_int = base_cl_int
        
        if ethnicity_distribution is None:
            self.ethnicity_distribution = {
                Ethnicity.EAST_ASIAN: 0.34,
                Ethnicity.EUROPEAN: 0.33,
                Ethnicity.AFRICAN: 0.33
            }
        else:
            self.ethnicity_distribution = ethnicity_distribution
        
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def generate(self) -> List[IndividualCharacteristics]:
        """가상 인구집단 생성
        
        Returns:
            IndividualCharacteristics 객체 리스트
        """
        population = []
        
        for i in range(self.n_subjects):
            # 1. 기본 인구통계학적 특성
            ethnicity = self._sample_ethnicity()
            gender = 'M' if np.random.random() < self.gender_ratio else 'F'
            age = np.random.randint(self.age_range[0], self.age_range[1] + 1)
            
            # 2. 신체 계측 (성별 고려)
            weight, height, bmi = self._sample_anthropometrics(gender, age)
            
            # 3. 유전자형 할당
            cyp2c19_genotype = self._sample_genotype('CYP2C19', ethnicity)
            cyp3a4_genotype = self._sample_genotype('CYP3A4', ethnicity)
            
            # 4. Activity Score 계산
            cyp2c19_as = self._calculate_activity_score('CYP2C19', cyp2c19_genotype)
            cyp3a4_as = self._calculate_activity_score('CYP3A4', cyp3a4_genotype)
            
            # 5. 종합 Activity Score (CYP2C19과 CYP3A4의 기하평균)
            combined_as = np.sqrt(cyp2c19_as * cyp3a4_as)
            
            # 6. 대사자 표현형 분류
            metabolizer_status = self._classify_metabolizer(combined_as)
            
            # 7. 생리학적 파라미터 생성
            phys_params = self._generate_physiological_params(
                weight, age, gender, combined_as
            )
            
            individual = IndividualCharacteristics(
                subject_id=i + 1,
                age=age,
                gender=gender,
                ethnicity=ethnicity,
                weight=weight,
                height=height,
                bmi=bmi,
                cyp2c19_genotype=cyp2c19_genotype,
                cyp3a4_genotype=cyp3a4_genotype,
                cyp2c19_activity_score=cyp2c19_as,
                cyp3a4_activity_score=cyp3a4_as,
                combined_activity_score=combined_as,
                metabolizer_status=metabolizer_status,
                phys_params=phys_params
            )
            
            population.append(individual)
        
        return population
    
    def _sample_ethnicity(self) -> Ethnicity:
        """민족 샘플링"""
        ethnicities = list(self.ethnicity_distribution.keys())
        probs = [self.ethnicity_distribution[e] for e in ethnicities]
        return np.random.choice(ethnicities, p=probs)
    
    def _sample_anthropometrics(
        self, gender: str, age: int
    ) -> Tuple[float, float, float]:
        """신체 계측치 샘플링
        
        성별과 나이에 따른 체중, 키 분포를 고려합니다.
        """
        # 성별에 따른 체중 조정
        if gender == 'M':
            weight_adj = 1.1  # 남성은 평균 10% 더 무거움
            height_mean = 175
            height_sd = 7
        else:
            weight_adj = 0.9
            height_mean = 162
            height_sd = 6
        
        # 나이에 따른 체중 조정 (중년 이후 약간 증가)
        age_factor = 1.0 + 0.005 * max(0, age - 40)
        
        weight = np.random.normal(
            self.weight_mean * weight_adj * age_factor,
            self.weight_sd
        )
        weight = max(40, min(weight, 150))  # 현실적인 범위로 제한
        
        height = np.random.normal(height_mean, height_sd)
        height = max(140, min(height, 200))
        
        bmi = weight / (height / 100) ** 2
        
        return weight, height, bmi
    
    def _sample_genotype(
        self, gene: str, ethnicity: Ethnicity
    ) -> Tuple[str, str]:
        """Hardy-Weinberg 평형을 사용한 유전자형 샘플링
        
        개별 대립유전자를 독립적으로 샘플링하여 유전자형을 생성합니다.
        JSON 파일에서 로드된 대립유전자 빈도를 사용합니다.
        """
        # JSON에서 로드된 데이터 사용
        allele_freqs = get_allele_frequencies(ethnicity, gene)
        alleles = list(allele_freqs.keys())
        probs = [allele_freqs[a] for a in alleles]
        
        # 두 대립유전자 독립 샘플링 (Hardy-Weinberg)
        allele1 = np.random.choice(alleles, p=probs)
        allele2 = np.random.choice(alleles, p=probs)
        
        # 정렬하여 일관된 표현 (*1/*2 = *1/*2, not *2/*1)
        return tuple(sorted([allele1, allele2]))
    
    def _calculate_activity_score(
        self, gene: str, genotype: Tuple[str, str]
    ) -> float:
        """유전자형에서 Activity Score 계산
        
        JSON 파일에서 로드된 Activity Score 데이터를 사용합니다.
        """
        score1 = get_activity_score(gene, genotype[0])
        score2 = get_activity_score(gene, genotype[1])
        return score1 + score2
    
    def _classify_metabolizer(self, activity_score: float) -> MetabolizerStatus:
        """Activity Score 기반 대사자 표현형 분류"""
        if activity_score <= 0.25:
            return MetabolizerStatus.POOR
        elif activity_score <= 1.0:
            return MetabolizerStatus.INTERMEDIATE
        elif activity_score <= 2.0:
            return MetabolizerStatus.NORMAL
        else:
            return MetabolizerStatus.ULTRA_RAPID
    
    def _generate_physiological_params(
        self, weight: float, age: int, gender: str, activity_score: float
    ) -> PhysiologicalParameters:
        """생리학적 파라미터 생성
        
        체중과 나이에 따른 장기 용적, 혈류량을 스케일링합니다.
        """
        # 체중 기반 스케일링
        weight_ratio = weight / 70.0
        
        # 혈장 용적 (약 4.7% of body weight)
        v_plasma = 0.047 * weight
        
        # 간 용적 (체중에 따라 스케일링, 약 2.5% of body weight)
        v_liver = 0.025 * weight
        
        # 간 혈류량 (표준: 90 L/h, 체중에 따라 스케일링)
        # 나이에 따른 감소 고려 (60세 이상 10% 감소)
        age_factor = 1.0 - 0.1 * max(0, (age - 60) / 20)
        q_liver = 90 * weight_ratio ** 0.75 * age_factor
        
        # 내재적 청소율 (Activity Score로 조정)
        # 추가적인 개인 간 변이 (CV 30%)
        individual_variability = np.random.lognormal(0, 0.3)
        cl_int = self.base_cl_int * activity_score * individual_variability
        
        # 신장 청소율 (선택적, 기본 0)
        cl_renal = 0.0
        
        return PhysiologicalParameters(
            body_weight=weight,
            v_plasma=v_plasma,
            v_liver=v_liver,
            q_liver=q_liver,
            cl_int=cl_int,
            cl_renal=cl_renal,
            activity_score=activity_score
        )
    
    def get_population_summary(
        self, population: List[IndividualCharacteristics]
    ) -> Dict:
        """인구집단 요약 통계"""
        
        # 기본 통계
        weights = [ind.weight for ind in population]
        ages = [ind.age for ind in population]
        activity_scores = [ind.combined_activity_score for ind in population]
        
        # 성별 분포
        n_male = sum(1 for ind in population if ind.gender == 'M')
        n_female = len(population) - n_male
        
        # 민족 분포
        ethnicity_counts = {}
        for eth in Ethnicity:
            count = sum(1 for ind in population if ind.ethnicity == eth)
            ethnicity_counts[eth.value] = count
        
        # 대사자 표현형 분포
        metabolizer_counts = {}
        for status in MetabolizerStatus:
            count = sum(1 for ind in population if ind.metabolizer_status == status)
            metabolizer_counts[status.value] = count
        
        return {
            'n_subjects': len(population),
            'demographics': {
                'age': {'mean': np.mean(ages), 'sd': np.std(ages), 
                        'min': min(ages), 'max': max(ages)},
                'weight': {'mean': np.mean(weights), 'sd': np.std(weights),
                          'min': min(weights), 'max': max(weights)},
                'gender': {'male': n_male, 'female': n_female,
                          'male_ratio': n_male / len(population)}
            },
            'ethnicity_distribution': ethnicity_counts,
            'metabolizer_distribution': metabolizer_counts,
            'activity_score': {
                'mean': np.mean(activity_scores),
                'sd': np.std(activity_scores),
                'min': min(activity_scores),
                'max': max(activity_scores)
            }
        }


def calculate_safety_margin(
    cmax_distribution: np.ndarray,
    toxic_threshold: float
) -> Dict:
    """안전 마진 계산
    
    Args:
        cmax_distribution: 개인별 Cmax 값 배열
        toxic_threshold: 독성 임계값 (ng/mL)
        
    Returns:
        안전 마진 분석 결과
    """
    n_total = len(cmax_distribution)
    n_exceeding = np.sum(cmax_distribution > toxic_threshold)
    
    return {
        'toxic_threshold': toxic_threshold,
        'n_total': n_total,
        'n_exceeding_threshold': int(n_exceeding),
        'percentage_exceeding': (n_exceeding / n_total) * 100,
        'percentage_safe': ((n_total - n_exceeding) / n_total) * 100,
        'cmax_max': float(np.max(cmax_distribution)),
        'cmax_95th_percentile': float(np.percentile(cmax_distribution, 95)),
        'safety_ratio': toxic_threshold / np.percentile(cmax_distribution, 95)
    }


if __name__ == "__main__":
    # Test population generation
    print("=== Population Generator Test ===\n")
    
    generator = PopulationGenerator(
        n_subjects=100,
        ethnicity_distribution={
            Ethnicity.EAST_ASIAN: 0.5,
            Ethnicity.EUROPEAN: 0.3,
            Ethnicity.AFRICAN: 0.2
        },
        age_range=(18, 65),
        gender_ratio=0.5,
        weight_mean=70,
        weight_sd=15,
        random_seed=42
    )
    
    population = generator.generate()
    summary = generator.get_population_summary(population)
    
    print(f"Generated {summary['n_subjects']} subjects\n")
    
    print("Demographics:")
    print(f"  Age: {summary['demographics']['age']['mean']:.1f} ± "
          f"{summary['demographics']['age']['sd']:.1f} years")
    print(f"  Weight: {summary['demographics']['weight']['mean']:.1f} ± "
          f"{summary['demographics']['weight']['sd']:.1f} kg")
    print(f"  Gender (M:F): {summary['demographics']['gender']['male']}:"
          f"{summary['demographics']['gender']['female']}")
    
    print("\nEthnicity Distribution:")
    for eth, count in summary['ethnicity_distribution'].items():
        print(f"  {eth}: {count} ({count/summary['n_subjects']*100:.1f}%)")
    
    print("\nMetabolizer Phenotypes:")
    for status, count in summary['metabolizer_distribution'].items():
        print(f"  {status}: {count} ({count/summary['n_subjects']*100:.1f}%)")
    
    print(f"\nActivity Score: {summary['activity_score']['mean']:.2f} ± "
          f"{summary['activity_score']['sd']:.2f}")
