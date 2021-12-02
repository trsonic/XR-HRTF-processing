function [Hoct]=cmplxsmooth(H,Noct,fs)

% fs=44100;
freq=(0:fs/(length(H)-1):fs/2)';
Noct=2*Noct;
H=H(1:end/2,:);

Noct=2*Noct;
% octave center frequencies
f1=1;
i=0;
while f1 < fs/2
    f1=f1*10^(3/(10*Noct));
    i=i+1;
    fc(i,:)=f1;
end

% octave edge frequencies
for i=0:length(fc)-1
    i=i+1;
    f1=10^(3/(20*Noct))*fc(i);
    fe(i,:)=f1;
end

% find nearest frequency edges
for i=1:length(fe)
    fe_p=find(freq>fe(i),1,'first');
    fe_m=find(freq<fe(i),1,'last');
    fe_0=find(freq==fe(i));
    if isempty(fe_0)==0
        fe(i)=fe_0;
    else
        p=fe_p-fe(i);
        m=fe(i)-fe_m;
        if p<m
            fe(i)=fe_p;
        else
           fe(i)=fe_m;
        end
    end
end

for i=1:length(fe)-1
    H_i=H(fe(i):fe(i+1),:);
    Hoct(i,1:size(H,2))=mean(H_i);
end

Habs=abs(H);
for i=1:length(fe)-1
    absHoct_i=Habs(fe(i):fe(i+1),:);
    absHoct(i,1:size(H,2))=mean(absHoct_i);
end
fc=fc(2:end);
Hoct=Hoct.*absHoct./abs(Hoct);
Hoct=interp1(fc,Hoct,freq,'spline');
Hoct=[Hoct;Hoct(end:-1:1,:)];



