function [x_fast,x_slow] = f(x,filter_width)

lpf_kernel=gaussian_kernel_1d(filter_width);
x_slow=truncated_conv_symmetric(x,lpf_kernel);
x_fast=x./x_slow;
