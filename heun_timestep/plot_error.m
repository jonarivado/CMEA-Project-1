clc
close all
clear all

FolderName = 'build/';
t1 = load([FolderName,'time1.txt']);
error1 = load([FolderName,'error1.txt']);
t2 = load([FolderName,'time2.txt']);
error2 = load([FolderName,'error2.txt']);
t3 = load([FolderName,'time3.txt']);
error3 = load([FolderName,'error3.txt']);
t4 = load([FolderName,'time4.txt']);
error4 = load([FolderName,'error4.txt']);
t5 = load([FolderName,'time5.txt']);
error5 = load([FolderName,'error5.txt']);
t6 = load([FolderName,'time6.txt']);
error6 = load([FolderName,'error6.txt']);
t7 = load([FolderName,'time7.txt']);
error7 = load([FolderName,'error7.txt']);
t8 = load([FolderName,'time8.txt']);
error8 = load([FolderName,'error8.txt']);
hold on
loglog(t1, error1);
loglog(t2, error2);
loglog(t3, error3);
loglog(t4, error4);
loglog(t5, error5);
loglog(t6, error6);
loglog(t7, error7);
loglog(t8, error8);

%legend('u','truesol','error');
xlabel('t');
ylabel('error');
