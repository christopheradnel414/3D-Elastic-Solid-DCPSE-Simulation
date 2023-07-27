clear all
clc
close all

file = fopen("geom_param.txt");
header = cell2mat( textscan(file,'%f',3,'collectoutput',1));
fclose(file);
file = fopen("geom_data.txt");
data = cell2mat( textscan(file, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t','collectoutput',1));
fclose(file);

E = header(1);
v = header(2);
h = header(3);

x = data(:,1);
y = data(:,2);
z = data(:,3);
nx = data(:,4);
ny = data(:,5);
nz = data(:,6);
dispX = data(:,7);
dispY = data(:,8);
dispZ = data(:,9);
forceX = data(:,10);
forceY = data(:,11);
forceZ = data(:,12);
bx = data(:,13);
by = data(:,14);
bz = data(:,15);

datasize = length(x);

counter = 1;
counter1 = 1;
for i = 1:datasize
    if isnan(forceX(i)) == false && isnan(forceY(i)) == false && isnan(forceZ(i)) == false
        if abs(forceX(i)) >= 0.000001 || abs(forceY(i)) >= 0.000001 || abs(forceZ(i)) >= 0.000001
            xneumann(counter) = x(i);
            yneumann(counter) = y(i);
            zneumann(counter) = z(i);
            nxboundary(counter) = nx(i);
            nyboundary(counter) = ny(i);
            nzboundary(counter) = nz(i);
            counter = counter + 1;
        end
    elseif isnan(dispX(i)) == false || isnan(dispY(i)) == false || isnan(dispZ(i)) == false
        xdirichlet(counter1) = x(i);
        ydirichlet(counter1) = y(i);
        zdirichlet(counter1) = z(i);
        counter1 = counter1 + 1;
    end
end

counter2 = 1;
for i = 1:datasize
    if abs(forceX(i)) <= 0.000001 && abs(forceY(i)) <= 0.000001 && abs(forceZ(i)) <= 0.000001
        xfree(counter2) = x(i);
        yfree(counter2) = y(i);
        zfree(counter2) = z(i);
        nxfree(counter2) = nx(i);
        nyfree(counter2) = ny(i);
        nzfree(counter2) = nz(i);
        counter2 = counter2 + 1;
    end
end

scattersize = 20;
scale = 1;
figure()
scatter3(x,y,z,scattersize,'filled')
hold on
% scatter3(xneumann,yneumann,zneumann,scattersize,'filled')
% hold on 
% quiver3(xneumann,yneumann,zneumann,nxboundary,nyboundary,nzboundary,scale,'black')
% hold on
scatter3(xdirichlet,ydirichlet,zdirichlet,scattersize,'filled')
hold on 
scatter3(xfree,yfree,zfree,scattersize,'filled')
hold on 
quiver3(xfree,yfree,zfree,nxfree,nyfree,nzfree,scale,'black')
daspect([1 1 1])
pbaspect([1 1 1])