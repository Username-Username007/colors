function signal_strength = color_signal(color, wavelen)
color = lower(color);
d = dictionary('r', {[570, 50]}, 'g', {[540, 40]}, 'b', {[450, 20]});
arg = d(color);
arg = arg{1};
dist = makedist('normal', mu = arg(1), sigma = arg(2));
pdf = @(x) dist.pdf(x) / dist.pdf(arg(1));
signal_strength = pdf(wavelen);
end