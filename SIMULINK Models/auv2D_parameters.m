clear
clc

m=30.48;
Iz=3.45;
Xu=-8.8065;
Yv=-65.5457;
Nr=-6.7352;
Xud=-0.93;
Yvd=-35.5;
Nrd=-35.5;

a1 = Xu/(m-Xud);
a2 = Yv/(m-Yvd);
a3 = Nr/(Iz-Nrd);
b1 = (m-Yvd)/(m-Xud);
b2 = (Xud-m)/(m-Yvd);
b3 = (Yvd-Xud)/(Iz-Nrd);
