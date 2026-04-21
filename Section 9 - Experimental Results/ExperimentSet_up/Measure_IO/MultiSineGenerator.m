function signal = MultiSineGenerator(Ndata,fmin,fmax,fs)
fr = fs/Ndata;

if fmin < fr
    fmin = fr;
else
    fmin = fr*floor(fmin/fr);
end

if fmax > fs/2
    fmax = fs/2;
else
    fmax = fr*floor(fmax/fr);
end

Nsines = round((fmax-fmin)/fr)+1; %Prevents issues where Nsines is not an integer

U = zeros(Ndata,1);

U(1+(fmin/fr):(Nsines+(fmin/fr)))=exp(1i*2*pi*rand(Nsines,1));
u=2*real(ifft(U)); u = u/std(u);

signal = u;
end
