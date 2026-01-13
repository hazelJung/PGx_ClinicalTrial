/**
 * Virtual Population PBPK Simulator - Frontend JavaScript
 * 이 파일을 수정하여 앱의 동작을 변경할 수 있습니다!
 */

// ============================================
// 전역 상태
// ============================================
let populationData = null;
let simulationData = null;

// ============================================
// DOM 요소 참조
// ============================================
document.addEventListener('DOMContentLoaded', function () {
    // 슬라이더 값 표시 업데이트
    setupSliderDisplays();

    // 버튼 이벤트 연결
    document.getElementById('generate-population').addEventListener('click', generatePopulation);
    document.getElementById('run-simulation').addEventListener('click', runSimulation);
    document.getElementById('fetch-pubchem').addEventListener('click', fetchPubChem);
    document.getElementById('toxic-threshold').addEventListener('change', updateSafetyReport);

    console.log('🧬 PBPK Simulator 초기화 완료');
});

// ============================================
// 슬라이더 디스플레이 설정
// ============================================
function setupSliderDisplays() {
    // 인원수 슬라이더
    const nSubjectsSlider = document.getElementById('n-subjects');
    const nSubjectsDisplay = document.getElementById('n-subjects-display');
    nSubjectsSlider.addEventListener('input', function () {
        nSubjectsDisplay.textContent = this.value;
    });

    // 성비 슬라이더
    const genderSlider = document.getElementById('gender-ratio');
    const genderDisplay = document.getElementById('gender-display');
    genderSlider.addEventListener('input', function () {
        genderDisplay.textContent = this.value + '%';
    });
}

// ============================================
// Step 1: 인구집단 생성
// ============================================
async function generatePopulation() {
    const btn = document.getElementById('generate-population');
    const statusDiv = document.getElementById('population-status');

    // 로딩 상태 표시
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> 생성 중...';
    statusDiv.innerHTML = '<span class="text-muted">인구집단을 생성하고 있습니다...</span>';

    // 파라미터 수집
    const params = {
        n_subjects: parseInt(document.getElementById('n-subjects').value),
        eth_asian: parseInt(document.getElementById('eth-asian').value),
        eth_european: parseInt(document.getElementById('eth-european').value),
        eth_african: parseInt(document.getElementById('eth-african').value),
        age_min: parseInt(document.getElementById('age-min').value),
        age_max: parseInt(document.getElementById('age-max').value),
        gender_ratio: parseInt(document.getElementById('gender-ratio').value),
        weight_mean: parseFloat(document.getElementById('weight-mean').value),
        weight_sd: parseFloat(document.getElementById('weight-sd').value),
        base_clint: parseFloat(document.getElementById('base-clint').value)
    };

    try {
        const response = await fetch('/api/generate-population', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params)
        });

        const data = await response.json();

        if (data.success) {
            populationData = data;

            // 디지털 트윈 그리드 업데이트
            renderDigitalTwins(data.individuals);

            // 인구집단 요약 업데이트
            renderPopulationSummary(data.summary);

            // 상태 메시지 & 시뮬레이션 버튼 활성화
            statusDiv.innerHTML = `<span class="status-success">✓ ${params.n_subjects}명 인구집단 생성 완료!</span>`;
            document.getElementById('run-simulation').disabled = false;
        } else {
            statusDiv.innerHTML = `<span class="status-error">❌ 오류: ${data.error}</span>`;
        }
    } catch (error) {
        statusDiv.innerHTML = `<span class="status-error">❌ 서버 오류: ${error.message}</span>`;
    } finally {
        btn.disabled = false;
        btn.innerHTML = '👥 인구집단 생성';
    }
}

// ============================================
// 디지털 트윈 그리드 렌더링
// ============================================
function renderDigitalTwins(individuals) {
    const grid = document.getElementById('digital-twin-grid');
    grid.innerHTML = '';

    // 최대 200명까지 표시
    const displayLimit = Math.min(200, individuals.length);

    for (let i = 0; i < displayLimit; i++) {
        const ind = individuals[i];
        const icon = document.createElement('div');

        // 성별 아이콘
        icon.textContent = ind.gender === 'M' ? '👨' : '👩';

        // 대사자 표현형에 따른 클래스
        const phenotypeClass = ind.metabolizer.toLowerCase();
        icon.className = `person-icon ${phenotypeClass}`;

        // 툴팁
        icon.title = `ID: ${ind.id} | ${ind.age}세 | ${ind.weight}kg | ${ind.ethnicity} | ${ind.metabolizer}`;

        grid.appendChild(icon);
    }

    // 더 많은 개인이 있으면 표시
    if (individuals.length > displayLimit) {
        const more = document.createElement('span');
        more.className = 'text-muted ms-2';
        more.textContent = `+${individuals.length - displayLimit}`;
        grid.appendChild(more);
    }
}

// ============================================
// 인구집단 요약 렌더링
// ============================================
function renderPopulationSummary(summary) {
    const container = document.getElementById('population-summary');
    const pieContainer = document.getElementById('metabolizer-pie');

    // 텍스트 요약
    const demo = summary.demographics;
    container.innerHTML = `
        <div class="stat-box">
            <div class="label">총 인원</div>
            <div class="value">${summary.n_subjects}명</div>
        </div>
        <div class="stat-box">
            <div class="label">나이</div>
            <div class="value">${demo.age.mean.toFixed(1)} ± ${demo.age.sd.toFixed(1)}세</div>
        </div>
        <div class="stat-box">
            <div class="label">체중</div>
            <div class="value">${demo.weight.mean.toFixed(1)} ± ${demo.weight.sd.toFixed(1)}kg</div>
        </div>
    `;

    // 파이 차트
    const metDist = summary.metabolizer_distribution;
    const pieData = [{
        values: [
            metDist['Poor Metabolizer (PM)'] || 0,
            metDist['Intermediate Metabolizer (IM)'] || 0,
            metDist['Normal Metabolizer (NM)'] || 0,
            metDist['Ultra-rapid Metabolizer (UM)'] || 0
        ],
        labels: ['PM', 'IM', 'NM', 'UM'],
        type: 'pie',
        hole: 0.4,
        marker: {
            colors: ['#E74C3C', '#F39C12', '#2ECC71', '#3498DB']
        },
        textinfo: 'label+percent'
    }];

    Plotly.newPlot(pieContainer, pieData, {
        margin: { l: 20, r: 20, t: 20, b: 20 },
        showlegend: false
    }, { responsive: true });
}

// ============================================
// Step 2: 시뮬레이션 실행
// ============================================
async function runSimulation() {
    if (!populationData) {
        alert('먼저 인구집단을 생성하세요!');
        return;
    }

    const btn = document.getElementById('run-simulation');
    const statusDiv = document.getElementById('simulation-status');

    // 로딩 상태
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> 시뮬레이션 중...';
    statusDiv.innerHTML = '<span class="text-muted">PBPK 시뮬레이션 실행 중...</span>';

    // 약물 파라미터
    const params = {
        drug_name: document.getElementById('drug-name').value,
        log_p: parseFloat(document.getElementById('log-p').value),
        f_u: parseFloat(document.getElementById('f-u').value),
        v_d: parseFloat(document.getElementById('v-d').value),
        k_a: parseFloat(document.getElementById('k-a').value),
        dose: parseFloat(document.getElementById('dose').value),
        bioavail: parseFloat(document.getElementById('bioavail').value),
        population: populationData.individuals
    };

    try {
        const response = await fetch('/api/run-simulation', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params)
        });

        const data = await response.json();

        if (data.success) {
            simulationData = data;

            // PK 곡선 그리기
            renderPKCurves(data);

            // 안전 마진 분석
            updateSafetyReport();

            statusDiv.innerHTML = `<span class="status-success">✓ 시뮬레이션 완료!</span>`;
        } else {
            statusDiv.innerHTML = `<span class="status-error">❌ 오류: ${data.error}</span>`;
        }
    } catch (error) {
        statusDiv.innerHTML = `<span class="status-error">❌ 서버 오류: ${error.message}</span>`;
    } finally {
        btn.disabled = false;
        btn.innerHTML = '💊 시뮬레이션 실행';
    }
}

// ============================================
// PK 곡선 렌더링
// ============================================
function renderPKCurves(data) {
    const container = document.getElementById('pk-curves-chart');

    const traces = [];

    // 개별 곡선 (일부만)
    if (data.individual_curves) {
        for (let i = 0; i < Math.min(30, data.individual_curves.length); i++) {
            traces.push({
                x: data.time,
                y: data.individual_curves[i],
                mode: 'lines',
                opacity: 0.3,
                line: { color: '#BDBDBD', width: 1 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }
    }

    // 90% CI
    traces.push({
        x: data.time.concat([...data.time].reverse()),
        y: data.ci_upper.concat([...data.ci_lower].reverse()),
        fill: 'toself',
        fillcolor: 'rgba(30, 58, 95, 0.2)',
        line: { color: 'transparent' },
        showlegend: true,
        name: '90% CI'
    });

    // 평균 곡선
    traces.push({
        x: data.time,
        y: data.mean_concentration,
        mode: 'lines',
        line: { color: '#1E3A5F', width: 3 },
        showlegend: true,
        name: '평균'
    });

    const layout = {
        xaxis: { title: '시간 (h)', gridcolor: '#E5E8EB' },
        yaxis: { title: '혈장 농도 (ng/mL)', gridcolor: '#E5E8EB' },
        margin: { l: 60, r: 30, t: 30, b: 50 },
        legend: { x: 0.02, y: 0.98 },
        plot_bgcolor: 'white'
    };

    Plotly.newPlot(container, traces, layout, { responsive: true });
}

// ============================================
// 안전 마진 분석
// ============================================
function updateSafetyReport() {
    if (!simulationData) return;

    const threshold = parseFloat(document.getElementById('toxic-threshold').value);
    const cmax = simulationData.cmax_distribution;

    // 임계값 초과 비율 계산
    const nExceeding = cmax.filter(c => c > threshold).length;
    const pctExceeding = (nExceeding / cmax.length * 100).toFixed(1);
    const cmax95 = percentile(cmax, 95).toFixed(1);
    const safetyRatio = (threshold / percentile(cmax, 95)).toFixed(2);

    // 보고서 렌더링
    const reportDiv = document.getElementById('safety-report');
    const alertClass = pctExceeding > 10 ? 'danger' : pctExceeding > 5 ? 'warning' : 'safe';

    reportDiv.innerHTML = `
        <div class="safety-alert ${alertClass}">
            <strong>안전 마진 분석 결과</strong><br>
            <span>${pctExceeding}% 인구가 독성 임계값 초과</span><br>
            <span>95th percentile Cmax: ${cmax95} ng/mL</span><br>
            <span>안전 비율: ${safetyRatio}</span>
        </div>
    `;

    // 히스토그램
    const histContainer = document.getElementById('cmax-histogram');
    Plotly.newPlot(histContainer, [{
        x: cmax,
        type: 'histogram',
        marker: { color: '#1E3A5F' },
        nbinsx: 30
    }], {
        xaxis: { title: 'Cmax (ng/mL)' },
        yaxis: { title: '빈도' },
        margin: { l: 50, r: 20, t: 20, b: 40 },
        shapes: [{
            type: 'line',
            x0: threshold, x1: threshold,
            y0: 0, y1: 1, yref: 'paper',
            line: { color: '#E74C3C', width: 2, dash: 'dash' }
        }]
    }, { responsive: true });
}

// ============================================
// PubChem 검색
// ============================================
async function fetchPubChem() {
    const drugName = document.getElementById('drug-name').value;
    const statusDiv = document.getElementById('pubchem-status');

    if (!drugName) {
        statusDiv.innerHTML = '<span class="text-warning">약물명을 입력하세요</span>';
        return;
    }

    statusDiv.innerHTML = '<span class="text-muted">검색 중...</span>';

    try {
        const response = await fetch(`/api/fetch-pubchem?drug_name=${encodeURIComponent(drugName)}`);
        const data = await response.json();

        if (data.found) {
            document.getElementById('log-p').value = data.log_p.toFixed(2);
            statusDiv.innerHTML = `<span class="status-success">✓ ${drugName} 발견 (MW: ${data.mw.toFixed(1)})</span>`;
        } else {
            statusDiv.innerHTML = '<span class="status-error">약물을 찾을 수 없습니다</span>';
        }
    } catch (error) {
        statusDiv.innerHTML = `<span class="status-error">오류: ${error.message}</span>`;
    }
}

// ============================================
// 유틸리티 함수
// ============================================
function percentile(arr, p) {
    const sorted = [...arr].sort((a, b) => a - b);
    const index = Math.ceil((p / 100) * sorted.length) - 1;
    return sorted[Math.max(0, index)];
}
