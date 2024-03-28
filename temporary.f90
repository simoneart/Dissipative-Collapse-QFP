!TERM IN W
(6*D2+12*D3-3*f)*f_of_t(t-1,j1,j2,j3,j4,j5,j6)

!TERM IN dq
-0.5/m/d*((Grd(j4)*(f_of_t(t-1,j1+1,j2,j3,j4,j5,j6)-f_of_t(t-1,j1-1,j2,j3,j4,j5,j6)) + &
Grd(j5)*(f_of_t(t-1,j1,j2+1,j3,j4,j5,j6)-f_of_t(t-1,j1,j2-1,j3,j4,j5,j6)) + &
Grd(j6)*(f_of_t(t-1,j1,j2,j3+1,j4,j5,j6)-f_of_t(t-1,j1,j2,j3-1,j4,j5,j6))))

!TERM IN d^2q
D4/d**2*((f_of_t(t-1,j1+1,j2,j3,j4,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1-1,j2,j3,j4,j5,j6)) + &
(f_of_t(t-1,j1,j2+1,j3,j4,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2-1,j3,j4,j5,j6)) + &
(f_of_t(t-1,j1,j2,j3+1,j4,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2,j3-1,j4,j5,j6)))

!TERM IN dp
((m*w**2*Grd(j1) + (4*D2+8*D3-f)*Grd(j4))*(f_of_t(t-1,j1,j2,j3,j4+1,j5,j6)-f_of_t(t-1,j1,j2,j3,j4-1,j5,j6)) + &
(m*w**2*Grd(j2) + (4*D2+8*D3-f)*Grd(j5))*(f_of_t(t-1,j1,j2,j3,j4,j5+1,j6)-f_of_t(t-1,j1,j2,j3,j4,j5-1,j6)) + &
(m*w**2*Grd(j3) + (4*D2+8*D3-f)*Grd(j6))*(f_of_t(t-1,j1,j2,j3,j4,j5,j6+1)-f_of_t(t-1,j1,j2,j3,j4,j5,j6-1)))/2.0/d

!TERM IN d^2p
((D1+D2*(Grd(j4)**2+Grd(j5)**2+Grd(j6)**2)+D3*Grd(j4)**2)*(f_of_t(t-1,j1,j2,j3,j4+1,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2,j3,j4-1,j5,j6)) + &
(D1+D2*(Grd(j4)**2+Grd(j5)**2+Grd(j6)**2)+D3*Grd(j5)**2)*(f_of_t(t-1,j1,j2,j3,j4,j5+1,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2,j3,j4,j5-1,j6)) + &
(D1+D2*(Grd(j4)**2+Grd(j5)**2+Grd(j6)**2)+D3*Grd(j6)**2)*(f_of_t(t-1,j1,j2,j3,j4,j5,j6+1) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2,j3,j4,j5,j6-1)) + &
2*D3*Grd(j4)*Grd(j5)*(f_of_t(t-1,j1,j2,j3,j4+1,j5+1,j6) - f_of_t(t-1,j1,j2,j3,j4-1,j5+1,j6) - f_of_t(t-1,j1,j2,j3,j4+1,j5-1,j6) + f_of_t(t-1,j1,j2,j3,j4-1,j5-1,j6))/4 + &
2*D3*Grd(j4)*Grd(j6)*(f_of_t(t-1,j1,j2,j3,j4+1,j5,j6+1) - f_of_t(t-1,j1,j2,j3,j4+1,j5,j6-1) - f_of_t(t-1,j1,j2,j3,j4-1,j5,j6+1) + f_of_t(t-1,j1,j2,j3,j4-1,j5,j6-1))/4 + &
2*D3*Grd(j5)*Grd(j6)*(f_of_t(t-1,j1,j2,j3,j4,j5+1,j6+1) - f_of_t(t-1,j1,j2,j3,j4,j5+1,j6-1) - f_of_t(t-1,j1,j2,j3,j4,j5-1,j6+1) + f_of_t(t-1,j1,j2,j3,j4,j5-1,j6-1))/4)/d**2 + &

!TERM OF IV ORDER
(3*D2*(f_of_t(t-1,j1+2,j2,j3,j4+2,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4+2,j5,j6) + f_of_t(t-1,j1-2,j2,j3,j4+2,j5,j6) -2*f_of_t(t-1,j1+2,j2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1-2,j2,j3,j4,j5,j6) + f_of_t(t-1,j1+1,j2,j3,j4-2,j5,j6) -2*f_of_t(t-1,j1,j2,j3,j4-2,j5,j6) + f_of_t(t-1,j1-2,j2,j3,j4-2,j5,j6)) + &
3*D2*(f_of_t(t-1,j1,j2+2,j3,j4,j5+2,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5+2,j6) + f_of_t(t-1,j1,j2-2,j3,j4,j5+2,j6) -2*f_of_t(t-1,j1,j2+2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2-2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2+2,j3,j4,j5-2,j6) -2*f_of_t(t-1,j1,j2,j3,j4,j5-2,j6) + f_of_t(t-1,j1,j2-2,j3,j4,j5-2,j6)) + &
3*D2*(f_of_t(t-1,j1,j2,j3+2,j4,j5,j6+2) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6+2) + f_of_t(t-1,j1,j2,j3-2,j4,j5,j6+2) -2*f_of_t(t-1,j1,j2,j3+2,j4,j5,j6) + &
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2,j3-2,j4,j5,j6) + f_of_t(t-1,j1,j2,j3+2,j4,j5,j6-2) -2*f_of_t(t-1,j1,j2,j3,j4,j5,j6-2) + f_of_t(t-1,j1,j2,j3-2,j4,j5,j6-2)))/16.0/d**4 + &
(D2*(f_of_t(t-1,j1+2,j2,j3,j4,j5+2,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5+2,j6) + f_of_t(t-1,j1-2,j2,j3,j4,j5+2,j6) -2*f_of_t(t-1,j1+2,j2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1-2,j2,j3,j4,j5,j6) + f_of_t(t-1,j1+2,j2,j3,j4,j5-2,j6) -2*f_of_t(t-1,j1,j2,j3,j4,j5-2,j6) + f_of_t(t-1,j1-2,j2,j3,j4,j5-2,j6)) + &
D2*(f_of_t(t-1,j1+2,j2,j3,j4,j5,j6+2) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6+2) + f_of_t(t-1,j1-2,j2,j3,j4,j5,j6+2) -2*f_of_t(t-1,j1+2,j2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1-2,j2,j3,j4,j5,j6) + f_of_t(t-1,j1+2,j2,j3,j4,j5,j6-2) -2*f_of_t(t-1,j1,j2,j3,j4,j5,j6-2) + f_of_t(t-1,j1-2,j2,j3,j4,j5,j6-2)) + &
D2*(f_of_t(t-1,j1,j2+2,j3,j4,j5,j6+2) - 2*f_of_t(t-1,j1,j2,j3,j4,j5,j6+2) + f_of_t(t-1,j1,j2-2,j3,j4,j5,j6+2) -2*f_of_t(t-1,j1,j2+2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2-2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2+2,j3,j4,j5,j6-2) -2*f_of_t(t-1,j1,j2,j3,j4,j5,j6-2) + f_of_t(t-1,j1,j2-2,j3,j4,j5,j6-2)) + &
D2*(f_of_t(t-1,j1,j2+2,j3,j4+2,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4+2,j5,j6) + f_of_t(t-1,j1,j2-2,j3,j4+2,j5,j6) -2*f_of_t(t-1,j1,j2+2,j3,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2-2,j3,j4,j5,j6) + f_of_t(t-1,j1,j2+2,j3,j4-2,j5,j6) -2*f_of_t(t-1,j1,j2,j3,j4-2,j5,j6) + f_of_t(t-1,j1,j2-2,j3,j4-2,j5,j6)) + &
D2*(f_of_t(t-1,j1,j2,j3+2,j4+2,j5,j6) - 2*f_of_t(t-1,j1,j2,j3,j4+2,j5,j6) + f_of_t(t-1,j1,j2,j3-2,j4+2,j5,j6) -2*f_of_t(t-1,j1,j2,j3+2,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2,j3-2,j4,j5,j6) + f_of_t(t-1,j1,j2,j3+2,j4-2,j5,j6) -2*f_of_t(t-1,j1,j2,j3,j4-2,j5,j6) + f_of_t(t-1,j1,j2,j3-2,j4-2,j5,j6)) + &
D2*(f_of_t(t-1,j1,j2,j3+2,j4,j5+2,j6) - 2*f_of_t(t-1,j1,j2,j3,j4,j5+2,j6) + f_of_t(t-1,j1,j2,j3-2,j4,j5+2,j6) -2*f_of_t(t-1,j1,j2,j3+2,j4,j5,j6) +&
4*f_of_t(t-1,j1,j2,j3,j4,j5,j6) -2*f_of_t(t-1,j1,j2,j3-2,j4,j5,j6) + f_of_t(t-1,j1,j2,j3+2,j4,j5-2,j6) -2*f_of_t(t-1,j1,j2,j3,j4,j5-2,j6) + f_of_t(t-1,j1,j2,j3-2,j4,j5-2,j6)))/16.0/d**4 + &
