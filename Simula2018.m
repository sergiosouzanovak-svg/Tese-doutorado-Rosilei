function [estM, estH, estC, estMD, estHD, estCD, pMc, pvM, pHc, pCOV, pC1c, pvF1, pCc, pvFT, i1] = Simula2018(R_phi_cell, R_tetra_cell, nf)
% SIMULA2018 Realiza a Analise Fatorial para as amostras geradas por Matriz1
%            comparando as matrizes de correlacao Phi e Tetracorica.
%
% Entradas:
%   R_phi_cell   = celula (1 x ns) com matrizes de correlacao Phi
%   R_tetra_cell = celula (1 x ns) com matrizes de correlacao Tetracorica
%   nf           = numero de fatores a extrair
%
% Saidas:
%   estM, estMD   = estatisticas do MSA [media, min, max, dp] (Phi x Tetracorica)
%   estH, estHD   = estatisticas das comunalidades [media, min, max, dp]
%   estC, estCD   = estatisticas das variancias explicadas [media, min, max, dp]
%   pMc, pvM      = p-valores dos testes t e F para o MSA
%   pHc, pCOV     = p-valores dos testes T2 de Hotelling e Box M (comunalidades)
%   pC1c, pvF1    = p-valores dos testes t e F para a explicacao do Fator 1
%   pCc, pvFT     = p-valores dos testes t e F para a explicacao Total
%   i1            = numero de amostras processadas

ns = length(R_phi_cell);
p  = size(R_phi_cell{1}, 1);

M   = zeros(1, ns);
MD  = zeros(1, ns);
H   = zeros(p, ns);
HD  = zeros(p, ns);
C   = zeros(nf, ns);
CD  = zeros(nf, ns);

i1 = 0;
for s = 1:ns
    R_phi   = R_phi_cell{s};
    R_tetra = R_tetra_cell{s};
    
    % Analise Fatorial para matriz Phi
    [MS_p, L1_p, h2_p, c3_p] = fator_matriz(R_phi, nf);
    
    % Analise Fatorial para matriz Tetracorica
    [MS_t, L1_t, h2_t, c3_t] = fator_matriz(R_tetra, nf);
    
    % Armazena resultados
    i1 = i1 + 1;
    M(i1)    = MS_p;
    MD(i1)   = MS_t;
    H(:, i1) = h2_p;
    HD(:, i1)= h2_t;
    C(:, i1) = c3_p(1:nf)';
    CD(:, i1)= c3_t(1:nf)';
end

% Estatisticas descritivas: [media, min, max, desvio-padrao]
estM  = [mean(M), min(M), max(M), std(M)];
estMD = [mean(MD), min(MD), max(MD), std(MD)];

estH  = [mean(H, 2), min(H, [], 2), max(H, [], 2), std(H, 0, 2)];
estHD = [mean(HD, 2), min(HD, [], 2), max(HD, [], 2), std(HD, 0, 2)];

estC  = [mean(C, 2), min(C, [], 2), max(C, [], 2), std(C, 0, 2)];
estCD = [mean(CD, 2), min(CD, [], 2), max(CD, [], 2), std(CD, 0, 2)];

% 1. Teste de variancias e medias para o MSA
VM  = sort([estM(4)^2, estMD(4)^2]);
FM  = VM(2) / VM(1);
pvM = 2 * (1 - fcdf(FM, ns - 1, ns - 1)); % bilateral
if pvM < 0.05
    tM  = (estM(1) - estMD(1)) / sqrt((estM(4)^2 / ns) + (estMD(4)^2 / ns));
    W1M = estM(4)^2 / ns;
    W2M = estMD(4)^2 / ns;
    GLM = round((W1M + W2M)^2 / ((W1M^2 / (ns - 1)) + (W2M^2 / (ns - 1))));
else
    Sp2 = ((ns - 1) * estM(4)^2 + (ns - 1) * estMD(4)^2) / (2 * ns - 2);
    tM  = (estM(1) - estMD(1)) / sqrt(Sp2 * (2 / ns));
    GLM = 2 * ns - 2;
end
pMc = 2 * (1 - tcdf(abs(tM), GLM)); % bilateral

% 2. Teste de variancias e medias para explicacao do Fator 1
VF1  = sort([estC(1, 4)^2, estCD(1, 4)^2]);
FF1  = VF1(2) / VF1(1);
pvF1 = 2 * (1 - fcdf(FF1, ns - 1, ns - 1));
if pvF1 < 0.05
    tC1  = (estC(1, 1) - estCD(1, 1)) / sqrt((estC(1, 4)^2 / ns) + (estCD(1, 4)^2 / ns));
    W1C1 = estC(1, 4)^2 / ns;
    W2C1 = estCD(1, 4)^2 / ns;
    GLC1 = round((W1C1 + W2C1)^2 / ((W1C1^2 / (ns - 1)) + (W2C1^2 / (ns - 1))));
else
    Sp2  = ((ns - 1) * estC(1, 4)^2 + (ns - 1) * estCD(1, 4)^2) / (2 * ns - 2);
    tC1  = (estC(1, 1) - estCD(1, 1)) / sqrt(Sp2 * (2 / ns));
    GLC1 = 2 * ns - 2;
end
pC1c = 2 * (1 - tcdf(abs(tC1), GLC1));

% 3. Teste de variancias e medias para explicacao Total
VFT  = sort([estC(nf, 4)^2, estCD(nf, 4)^2]);
FFT  = VFT(2) / VFT(1);
pvFT = 2 * (1 - fcdf(FFT, ns - 1, ns - 1));
if pvFT < 0.05
    tC  = (estC(nf, 1) - estCD(nf, 1)) / sqrt((estC(nf, 4)^2 / ns) + (estCD(nf, 4)^2 / ns));
    W1C = estC(nf, 4)^2 / ns;
    W2C = estCD(nf, 4)^2 / ns;
    GLC = round((W1C + W2C)^2 / ((W1C^2 / (ns - 1)) + (W2C^2 / (ns - 1))));
else
    Sp2 = ((ns - 1) * estC(nf, 4)^2 + (ns - 1) * estCD(nf, 4)^2) / (2 * ns - 2);
    tC  = (estC(nf, 1) - estCD(nf, 1)) / sqrt(Sp2 * (2 / ns));
    GLC = 2 * ns - 2;
end
pCc = 2 * (1 - tcdf(abs(tC), GLC));

% 4. Teste de igualdade de matrizes de covariancia das comunalidades (Box M)
SH   = cov(H');
SHD  = cov(HD');
S_pool = ((ns - 1) * SH + (ns - 1) * SHD) / (2 * (ns - 1));
MCOV = 2 * (ns - 1) * log(det(S_pool)) - (ns - 1) * (log(det(SH)) + log(det(SHD)));
CI   = 1 - (((2 * p^2 + 3 * p - 1) * 3) / (6 * (p + 1) * 3 * ns));
Q2   = MCOV * CI;
GL_box = 0.5 * p * (p + 1);
pCOV = 1 - chi2cdf(real(Q2), GL_box);

% 5. Teste do vetor medio de comunalidades (T2 de Hotelling)
m1_h = mean(H, 2)';
m2_h = mean(HD, 2)';
diff_h = m1_h - m2_h;
if pCOV < 0.05
    T2  = diff_h * inv((SH / ns) + (SHD / ns)) * diff_h';
    pHc = 1 - chi2cdf(real(T2), p);
else
    Sp_h = ((ns - 1) * SH + (ns - 1) * SHD) / (2 * ns - 2);
    T2   = diff_h * inv((2 / ns) * Sp_h) * diff_h';
    F_hot = T2 * (2 * ns - p - 1) / ((2 * ns - 2) * p);
    pHc  = 1 - fcdf(real(F_hot), p, 2 * ns - p - 1);
end
