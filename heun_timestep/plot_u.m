clc
close all
clear all

FolderName = 'build/';
u = load([FolderName,'solution8.txt']);
t = load([FolderName,'time8.txt']);
error = load([FolderName,'error8.txt']);
truesol = load([FolderName,'truesol8.txt']);
hold on
plot(t, u);
plot(t,truesol);
plot(t,error);
legend('u','truesol','error');
xlabel('t');
ylabel('u(t)');
