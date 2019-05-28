tic

%% Limpeza

clc; clear all; close all;

%% Entradas

[y,Fs] = audioread('C:\Users\msnda\Desktop\110118\110118_001.WAV');

yTx = y(:,5);
y1 = y(:,1);
num_amostras = 256;
i1 = 1;
i2 = i1 + num_amostras - 1;
rot = 0;

while i2 < length(yTx)
    yblk = yTx(i1:i2);
    iTrg = find(yblk > 0.05); %0.0002
    while isempty(iTrg)
        i1 = i1 + num_amostras;
        i2 = i1 + num_amostras - 1;
        if i2 > length(yTx); break; end
        
        yblk = yTx(i1:i2);
        
        iTrg = find(yblk > 0.05); %0.0002

    end
    if isempty(iTrg); break; end
    i11 = i1 + iTrg(1) + 300;

    i1 = i1 + 24000;

    i2 = i1 + num_amostras - 1;
        if i2 > length(yTx); break; end

    
    sinal = y(i11:i11 + num_amostras - 1,1:4);
    rot = rot + 1;
    
    v_som = 1491.24; % velocidade do som no meio em m/s
    d = 4 * 10^-2; % distancia entre os hidrofones em m
    num_canais = 4; % numero de hidrofones
    fs = Fs; % frequencia de amostragem em Hz
    
    
    % montagem da matriz de atrasos
    d_theta = 1; % resolucao angular em graus
    az = length(0:d_theta:180); % numero de angulos
    d_phi = 1; % resolucao angular em graus
    el = length(0:d_theta:180); % numero de angulos
    delta_m = zeros(num_canais,az);
    
    % matriz de coordenadas arranjo
        coord = [0 0 d
                 0 0 0
                 0 0 -d
                -d 0 0];
        
        for n = 1:az
            for m = 1:el
                for k = 1:num_canais
                    theta = (n - 1) * d_theta;
                    phi = (m - 1) * d_phi;
                    s = [cosd(theta)*sind(phi) sind(theta)*sind(phi) -cosd(phi)];
                    delta_m(k,(n-1)*el + m) = (dot(s,coord(k,:)))/v_som;
                    par((n-1)*el+m,:) = [theta phi]; % pares de angulos para consulta
                end
            end
        end

      
        for k = 0 : num_amostras/2 - 1
            delay(k + 1,:,:) = exp((1j * 2 * pi * fs)/(num_amostras) * k * delta_m);
        end
        
        FFT_sinal = fft(sinal,num_amostras);
        
        for a = 1 : az*el
            aux = real(ifft(FFT_sinal(1:end/2,:) .* delay(:,:,a)));
            soma = (1/num_canais)*sum(aux,2);
            dados(rot,a) = (((1/num_amostras/2)*(sum(soma.^2)))^0.5);
        end
end


dados = reshape(dados,[rot,az,el]);
dados = dados(1,:,:);
dados = reshape(dados,[az,el]);

dados = dados/max(max(dados));


imagesc(0:az-1,0:el-1,dados)
set(gca,'YDir','normal')
xlabel('Azimute (º)')
ylabel('Elevação (º)')


toc