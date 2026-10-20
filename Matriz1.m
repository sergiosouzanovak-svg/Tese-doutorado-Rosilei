function [amostras_normais, amostras_dicot, R_phi_cell, R_tetra_cell] = Matriz1(R, n, ns, per)
% MATRIZ1 Gera ns amostras multivariadas normalmente distribuidas pelo metodo
%         de Monte Carlo, dicotomiza sob tres condicoes de zc e calcula as
%         matrizes de correlacao Phi e Tetracorica para cada amostra valida.
%
% Entradas:
%   R   = matriz de correlacao teorica/populacional (ex: gerada por matrizc1)
%   n   = tamanho da amostra (ex: 5, 8, 10, 12, 15, 18, 20 vezes o numero de variaveis)
%   ns  = numero de simulacoes validas (padrao = 100)
%   per = ponto de dicotomizacao:
%         1 -> P(z <= zc) = 0.25 (zc = -0.6745)
%         2 -> P(z <= zc) = 0.50 (zc =  0.0000)
%         3 -> P(z <= zc) = 0.75 (zc =  0.6745)
%
% Saidas:
%   amostras_normais = celula (1 x ns) com as matrizes X (n x p) normais
%   amostras_dicot   = celula (1 x ns) com as matrizes XD (n x p) dicotomizadas (0/1)
%   R_phi_cell       = celula (1 x ns) com as matrizes de correlacao Phi (corrcoef de XD)
%   R_tetra_cell     = celula (1 x ns) com as matrizes de correlacao Tetracorica
%
% Compativel com MATLAB R12 em diante.

if nargin < 4
    error('Uso: Matriz1(R, n, ns, per)');
end

p = size(R, 1);
m = zeros(1, p);

% Ponto de corte zc conforme o percentil teorico da normal padrao
if per == 1
    ZC = -0.6745; % P(z <= zc) = 0.25 (25% de 0 e 75% de 1)
elseif per == 2
    ZC = 0.0000;  % P(z <= zc) = 0.50 (50% de 0 e 50% de 1)
elseif per == 3
    ZC = 0.6745;  % P(z <= zc) = 0.75 (75% de 0 e 25% de 1)
else
    error('Parametro per deve ser 1, 2 ou 3.');
end

% Pre-alocacao das celulas de saida
amostras_normais = cell(1, ns);
amostras_dicot   = cell(1, ns);
R_phi_cell       = cell(1, ns);
R_tetra_cell     = cell(1, ns);

i1 = 0; % Contador de simulacoes validas aceitas

while i1 < ns
    % 1. Geracao da amostra normal multivariada pelo metodo de Monte Carlo
    X = mvnrnd(m, R, n);
    
    % 2. Dicotomizacao em 0 e 1 de acordo com o limiar ZC
    XD = zeros(n, p);
    for ii = 1:n
        for jj = 1:p
            if X(ii, jj) < ZC
                XD(ii, jj) = 0;
            else
                XD(ii, jj) = 1;
            end
        end
    end
    
    % 3. Verificacao de variancia nula (colunas constantes)
    sd_cols = std(XD);
    if any(sd_cols == 0) || any(~isfinite(sd_cols))
        continue; % descarta e reamostra
    end
    
    % 4. Calculo das matrizes de correlacao Phi e Tetracorica
    R_phi   = corrcoef(XD);          % Phi equivale a Pearson sobre {0,1}
    R_tetra = tetracorica(XD);        % Matriz de correlacao tetracorica
    
    % Validar finitude das matrizes e autovalores
    if any(~isfinite(R_phi(:))) || any(~isfinite(R_tetra(:)))
        continue;
    end
    
    eig_phi   = eig(R_phi);
    eig_tetra = eig(R_tetra);
    if any(eig_phi <= 0) || any(eig_tetra <= 0)
        % Se a matriz nao for positiva definida, ajusta ou reamostra
        continue;
    end
    
    % 5. Amostra valida aceita
    i1 = i1 + 1;
    amostras_normais{i1} = X;
    amostras_dicot{i1}   = XD;
    R_phi_cell{i1}       = R_phi;
    R_tetra_cell{i1}     = R_tetra;
end
