function output = echo_gen(s, Fs, delay, amp)
echo_sample = delay*Fs;
echo_fmat = [zeros(echo_sample,1);s];
echo_smat = [s*amp;zeros(echo_sample,1)];
output = echo_fmat + echo_smat;
a = max(abs(output));
if a > 1
    output=output/a;
end
end