function [bits, estado_final] = LFSR_seq(m, taps, estado_ini, L)
% LFSR_seq  Genera una secuencia pseudoaleatoria usando un LFSR.
%
%   [bits, estado_final] = LFSR_seq(m, taps, estado_ini, L)
%
%   ENTRADAS:
%       m          -> longitud del registro (número de flip-flops)
%       taps       -> posiciones de realimentación (vector, ej: [5 4 2 1])
%       estado_ini -> vector fila con el estado inicial (1 x m), no todo ceros
%       L          -> longitud de la secuencia a generar
%
%   SALIDAS:
%       bits         -> secuencia generada (1 x L) con valores 0/1
%       estado_final -> estado del registro después de generar L bits
%
%   NOTAS:
%       - Se asume que el bit de salida se toma del último flip-flop (estado(end)).
%       - Las posiciones en 'taps' se refieren a índices de 'estado' empezando en 1.
%         (estado(1) es el bit más a la izquierda).
%
%   EJEMPLO:
%       m = 7; taps = [7 6 3 1]; estado_ini = [1 0 0 0 0 0 0]; L = 1000000;
%       [bits, est_fin] = LFSR_seq(m, taps, estado_ini, L);

    % --- Comprobaciones básicas ---
    if length(estado_ini) ~= m
        error('La longitud de estado_ini (%d) debe ser igual a m (%d).', ...
               length(estado_ini), m);
    end

    if all(estado_ini == 0)
        error('El estado inicial no puede ser todo ceros para un LFSR.');
    end

    % Si no se pasa L o viene vacío, generamos el período máximo teórico
    if nargin < 4 || isempty(L)
        L = 2^m - 1;
    end

    % --- Inicializaciones ---
    estado = estado_ini(:).';      % aseguramos que sea vector fila
    bits   = zeros(1, L);          % prealocación de la secuencia

    % --- Bucle principal del LFSR ---
    for n = 1:L
        % Salida del LFSR: último flip-flop
        salida   = estado(end);
        bits(n)  = salida;

        % Nuevo bit de entrada: XOR de los taps
        nuevo_bit = 0;
        for k = 1:length(taps)
            nuevo_bit = xor(nuevo_bit, estado(taps(k)));
        end

        % Desplazamiento a la derecha + inserción del nuevo_bit al inicio
        estado = [nuevo_bit, estado(1:end-1)];
    end

    % Estado final después de L iteraciones
    estado_final = estado;
end
