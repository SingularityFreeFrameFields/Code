function [Fout] = ComplexToLong(Fin, Gmag, narrowBand)

tempF = Fin(1:length(Fin) / 2) + 1i * Fin(length(Fin) / 2 + 1 : end);
Fout = tempF(Gmag(narrowBand) ~= 0);
end

