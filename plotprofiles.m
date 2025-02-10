function [] = plotprofiles(turb, can)

% Feb 2025
%HigherOrderClosure_v2025

nfig=1;

figure(nfig)
clf
plot(turb.uw(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.duw(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
xlabel('wu')
ylim([0 300])
legend('wu(z)/u*^2','dwu/u*^2')
title('momentum transport')

nfig=nfig+1;
figure(nfig);
clf
plot(turb.u(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)  % I think this is u/u*
xlabel('u/u*')  
ylim([0 300])
legend('u')

nfig=nfig+1;
figure(nfig);
clf
plot(turb.u(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)  % I think this is u/u*
hold on
plot(turb.du(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)  % I think this is u/u*
xlabel('u/u*')  
ylim([0 300])
legend('u','du')

nfig=nfig+1;
figure(nfig);
clf
plot(turb.wwu(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.wuu(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.wvv(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.w3(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)

xlabel(' xxx/u*^3 ')
ylim([0 300])
legend('wwu', 'wuu','wvv','www')
title('third order transport terms')


nfig=nfig+1;
figure(nfig);
clf
title ('turbulence variance profiles')
plot(turb.w2(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.u2(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.v2(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)
hold on
plot(turb.qq(1:can.sze3),1:can.sze3,'.-','MarkerSize', 10)

xlabel('xx/u*^2')
ylim([0 300])
legend('ww','uu','vv','qq, tke')


LAIcum=cumsum(can.lai_z);

nfig=nfig+1;
figure(nfig);
clf
plot(LAIcum(1:can.sze),1:can.sze,'.-','MarkerSize', 10)
xlabel('cumulative LAI')

nfig=nfig+1;
figure(nfig);
clf
plot(can.lai_z(1:can.sze),1:can.sze,'.-','MarkerSize', 10)
hold on
plot(can.pad(1:can.sze),1:can.sze,'.-','MarkerSize', 10)
hold on
xlabel('LAI Profiles')
ylabel('layers')
legend('lai = f(z)', 'leaf area density, m2 m-3')


cumLAI=can.lai-cumsum((can.lai_z(1:can.sze)));

nfig=nfig+1;
figure(nfig);
clf
title ('turbulence variance profiles')
plot(turb.w2(1:can.sze),cumLAI,'.-','MarkerSize', 10)
hold on
plot(turb.u2(1:can.sze),cumLAI,'.-','MarkerSize', 10)
hold on
plot(turb.v2(1:can.sze),cumLAI,'.-','MarkerSize', 10)
hold on
plot(turb.qq(1:can.sze),cumLAI,'.-','MarkerSize', 10)
ylabel('cumulative LAI')
xlabel('variance m2 s-2')
legend('ww','uu','vv','qq')
set(gca, 'YDir','reverse')

end