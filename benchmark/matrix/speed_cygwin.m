
%oodles of speed tests in MFLOPS of the various alrogithms
% for matrix multiplication as wells as a slew of data types
% these were done on a 933 MHz Pentium III under cygwin
% compiled under gcc-3.2.1

Infinity=Inf;
R=[2 5 20 24 31 32 36 48 64 73 96 97 127 128 129 163 191 192 229 255 256 257 319 320 321 417 479 480 511 512 ];
basic_mulD=[9.945 53.7621 82.28 81.4135 82.171 82.5955 81.5286 75.4975 68.1654 72.4279 68.3441 71.233 71.0936 66.052 68.6254 65.2467 36.7698 34.845 35.3335 35.3736 10.5917 35.3269 35.5161 33.7597 35.3567 34.1553 33.1593 33.036 32.9718 10.0835  ];

basic_mulC=[9.09827 36.1678 48.7982 51.7579 51.3569 51.1111 47.7757 48.3132 44.7392 47.6044 43.1579 45.5625 46.6204 42.5278 45.796 43.0385 28.3103 27.0536 27.5752 27.6529 9.51224 27.5897 27.5625 26.2144 27.4919 26.3296 25.9425 25.7378 25.4534 9.19937  ]*2;

basic_mulZ=[18.04 72.98 95.3251 92.9772 96.4616 96.5595 86.0824 82.3609 84.2018 86.0746 81.4721 83.0883 54.4873 38.0005 53.0047 36.643 39.2004 35.3563 39.1651 38.3274 16.6307 36.5143 39.1991 16.8869 37.8013 35.7674 35.7726 24.1 35.7777 16.177  ]*2;

at_mulD=[21.468 162.822 583.515 585.63 612.416 642.19 611.53 658.887 715.828 662.267 658.887 649.915 583.624 572.662 567.134 591.291 600.598 597.605 607.089 606.058 613.918 614.813 615.754 627.139 616.087 623.757 625.333 629.257 626.445 624.995 ];

at_mulC=[9.09827 63.5655 167.718 172.631 175.574 181.008 182.524 183.543 357.914 299.965 265.992 241.243 322.699 340.007 308.391 260.987 320.592 357.808 351.268 329.773 361.529 331.212 332.515 356.901 319.576 360.531 347.38 367.568 363.453 341.956 ]*2;

at_mulZ=[19.7031 113.384 322.639 325.654 354.979 343.487 359.554 329.444 324.002 318.428 299 302.061 313.816 312.316 305.817 314.071 316.947 321.951 323.747 326.526 327.161 302.443 291.627 322.44 319.576 333.004 334.94 341.004 284.695 234.135 ]*2;



clf reset;
colordef white;

subplot(2,2,1);
hold on;
semilogx(R, basic_mulD,'k-');
semilogx(R, matlab_mulD,'k*-');
semilogx(R, at_mulD,'k.-');
title('C=A*B--Double Matrix Multiplication 933 MHz Pentium III');
xlabel('NxN');
ylabel('MFLOPS');
legend('Basic','Matlab 5.3','ATLAS 3.4');
hold off;

subplot(2,2,2);
hold on;
semilogx(R, basic_mulC,'k-');
semilogx(R, at_mulC,'k.-');
title('C=A*B--complex<float> Matrix Multiplication 933 MHz Pentium III');
xlabel('NxN');
ylabel('MFLOPS');
legend('Basic','ATLAS 3.4');
hold off;

subplot(2,2,3);
hold on;
semilogx(R, basic_mulZ,'k-');
semilogx(R, matlab_mul*2,'k*-');
semilogx(R, at_mulZ,'k.-');
title('C=A*B--complex<double> Matrix Multiplication 933 MHz Pentium III');
xlabel('NxN');
ylabel('MFLOPS');
legend('Basic','Matlab 5.3','ATLAS 3.4');
hold off;
