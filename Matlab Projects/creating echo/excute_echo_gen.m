function excute_echo_gen()
load splat
output = echo_gen(y, Fs, 0.25, 0.6);
dt = 1/Fs; 
t = 0:dt:dt*(length(output)-1);
plot(t, output)
sound(y,Fs)
pause(2);
sound(output,Fs)
end