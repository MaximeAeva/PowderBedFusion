%% Estimation des paramètres de fusion.
clear all; close all;


duree=15000;
distancex=300; %nombre de cellule de la grille
distancey=600;
taille_faisceau=66;
vitesse_faisceau=30;
P_laser=10;%en W




T=zeros(distancex,distancey);%création de la matrice
Position_particules=zeros(distancex,distancey);
dx=10e-6;%en m
dV=dx^3;


capa_air=1/0.026;
diff_air=20;

for i=1:distancex
    for j=1:distancey
        capa(i,j)=capa_air;
        diff(i,j)=diff_air;
    end
end

lc=4*taille_faisceau;
kc=lc;

%***Position éléments***

%     Ajout de billes

capa_poudre=1/427;
diff_poudre=173;
T_fusion_poudre=962;

nbparticules=2000;

diam=8;          % diamètre sphère
rad=diam/2;     % rayon

icentre = randi([1+diam distancex-diam],1,nbparticules);
jcentre = randi([1+diam distancey-diam],1,nbparticules);

%icentre=[distancex/2+2*rad-2 distancex/2 distancex/2-2*rad+2 distancex/2-4*rad+4 distancex/2+4*rad-4];   % i-centre
%jcentre=[distancey/2 distancey/2 distancey/2 distancey/2 distancey/2];      % j-centre
for n=1:length(icentre)
for i=1:distancex
    for j=1:distancey
        dist2=(i+0.5-icentre(n))^2 + (j-jcentre(n))^2;
        if dist2 <= rad^2 
         capa(i,j)=capa_poudre;
         diff(i,j)=diff_poudre;
         Position_particules(i,j)=1;
        end
    end
end
end

%***__________***

figure('units','normalized','outerposition',[0 0 1 1]);
for tmp=0:duree
    
%Source
    if tmp<(distancey-2*(taille_faisceau+1))/vitesse_faisceau
        tmps=tmp+((taille_faisceau+1)/vitesse_faisceau);
        for l=1:lc
            for k=1:kc
                dist=(((lc/2)-l))^2+(((kc/2)-k))^2;
                if dist <= taille_faisceau^2 
                    T(((lc/2)-l)+distancex/2,floor(vitesse_faisceau*tmps-((kc/2)-k)))=capa(((lc/2)-l)+distancex/2,floor(vitesse_faisceau*tmps-((kc/2)-k)))*...
                        (T(((lc/2)-l)+distancex/2,floor(vitesse_faisceau*tmps-((kc/2)-k)))+...
                        (P_laser)*exp(-dist/(2*(taille_faisceau/5)^2)));
                end
            end
        end
    end
    
%Affichage
    LB=flipud(lbmap(256,'RedBlue')); 
    colormap(LB);
    %surf(T);
    surf(Position_particules);
    shading flat;
    %caxis([0 T_fusion_poudre]);
    axis([1 distancex 1 distancey]);
    axis image; axis xy
    axis off;
    timestep=int2str(tmp);
    title(['temperature à t = ',timestep]);
    colorbar
    pause(0.00001);
    
%Condition haut gauche
    T(1,distancey)=(4*capa(1,distancey)*T(1,distancey)+...
                    2*diff(2,distancey)*T(2,distancey)+2*diff(1,distancey-1)*T(1,distancey-1)+...
                    diff(2,distancey-1)*T(2,distancey-1))/...
                    (4*capa(1,distancey)+2*(diff(2,distancey)+diff(1,distancey-1))+diff(2,distancey-1));
                
%Condition haut droit
    T(distancex,distancey)=(4*capa(distancex,distancey)*T(distancex,distancey)+...
                            2*diff(distancex-1,distancey)*T(distancex-1,distancey)+2*diff(distancex,distancey-1)*T(distancex,distancey-1)+...
                            diff(distancex-1,distancey-1)*T(distancex-1,distancey-1))/...
                            (4*capa(distancex,distancey)+2*(diff(distancex-1,distancey)+diff(distancex,distancey-1))+diff(distancex-1,distancey-1));
                        
%Condition bas gauche
    T(1,1)=(4*capa(1,1)*T(1,1)+...
            2*diff(2,1)*T(2,1)+2*diff(1,2)*T(1,2)+...
            diff(2,2)*T(2,2))/...
            (4*capa(1,1)+2*(diff(2,1)+diff(1,2))+diff(2,2));
        
%Condition bas droit
    T(distancex,1)=(4*capa(distancex,1)*T(distancex,1)+...
                    2*diff(distancex-1,1)*T(distancex-1,1)+2*diff(distancex,2)*T(distancex,2)+...
                    diff(distancex-1,2)*T(distancex-1,2))/...
                    (4*capa(distancex,1)+2*(diff(distancex-1,1)+diff(distancex,2))+diff(distancex-1,2));
                
%Condition gauche
    T(1,2:distancey-1)=(4*capa(1,2:distancey-1).*T(1,2:distancey-1)+...
                        2*diff(1,1:distancey-2).*T(1,1:distancey-2)+2*diff(1,3:distancey).*T(1,3:distancey)+2*diff(2,2:distancey-1).*T(2,2:distancey-1)+...
                        diff(2,1:distancey-2).*T(2,1:distancey-2)+diff(2,3:distancey).*T(2,3:distancey))/...
                        (4*capa(1,2:distancey-1)+2*(diff(1,1:distancey-2)+diff(1,3:distancey)+diff(2,2:distancey-1))+diff(2,1:distancey-2)+diff(2,3:distancey));
                    
%Condition droite
    T(distancex,2:distancey-1)=(4*capa(distancex,2:distancey-1).*T(distancex,2:distancey-1)+...
                        2*diff(distancex,1:distancey-2).*T(distancex,1:distancey-2)+2*diff(distancex,3:distancey).*T(distancex,3:distancey)+2*diff(distancex-1,2:distancey-1).*T(distancex-1,2:distancey-1)+...
                        diff(distancex-1,1:distancey-2).*T(distancex-1,1:distancey-2)+diff(distancex-1,3:distancey).*T(distancex-1,3:distancey))/...
                        (4*capa(distancex,2:distancey-1)+2*(diff(distancex,1:distancey-2)+diff(distancex,3:distancey)+diff(distancex-1,2:distancey-1))+diff(distancex-1,1:distancey-2)+diff(distancex-1,3:distancey));
                    
%Condition haut
    T(2:distancex-1,distancey)=(4*capa(2:distancex-1,distancey)'.*T(2:distancex-1,distancey)+...
                        2*diff(1:distancex-2,distancey)'.*T(1:distancex-2,distancey)+2*diff(3:distancex,distancey)'.*T(3:distancex,distancey)+2*diff(2:distancex-1,distancey-1)'.*T(2:distancex-1,distancey-1)+...
                        diff(1:distancex-2,distancey-1)'.*T(1:distancex-2,distancey-1)+diff(3:distancex,distancey-1)'.*T(3:distancex,distancey-1))/...
                        (4*capa(2:distancex-1,distancey)'+2*(diff(1:distancex-2,distancey)'+diff(3:distancex,distancey)'+diff(2:distancex-1,distancey-1)')+diff(3:distancex,distancey-1)'+diff(1:distancex-2,distancey-1)');
                    
%Condition bas
    T(2:distancex-1,1)=(4*capa(2:distancex-1,1)'.*T(2:distancex-1,1)+...
                        2*diff(1:distancex-2,1)'.*T(1:distancex-2,1)+2*diff(3:distancex,1)'.*T(3:distancex,1)+2*diff(2:distancex-1,2)'.*T(2:distancex-1,2)+...
                        diff(1:distancex-2,2)'.*T(1:distancex-2,2)+diff(3:distancex,2)'.*T(3:distancex,2))/...
                        (4*capa(2:distancex-1,1)'+2*(diff(1:distancex-2,1)'+diff(3:distancex,1)'+diff(2:distancex-1,2)')+diff(1:distancex-2,2)'+diff(3:distancex,2)');
                    
%Itérations
    for i=2:distancex-1
        for j=2:distancey-1
            T(i,j)=(4*capa(i,j)*T(i,j)+...%Inertie de la cellule
            2*diff(i-1,j)*T(i-1,j)+2*diff(i+1,j)*T(i+1,j)+2*diff(i,j-1)*T(i,j-1)+2*diff(i,j+1)*T(i,j+1)+...%Chaleur transmise par les cellules en contact direct
            diff(i-1,j-1)*T(i-1,j-1)+diff(i-1,j+1)*T(i-1,j+1)+diff(i+1,j-1)*T(i+1,j-1)+diff(i+1,j+1)*T(i+1,j+1))/...%chaleur transmise par les cellules en contact indirect
            (4*capa(i,j)+2*(diff(i-1,j)+diff(i+1,j)+diff(i,j-1)+diff(i,j+1))+diff(i-1,j-1)+diff(i-1,j+1)+diff(i+1,j-1)+diff(i+1,j+1));
        end
    end
    %Force de Marangoni
    for i=2:distancex-1
        for j=2:distancey-1
            if capa(i,j)~=capa_air && T(i,j)>10*T_fusion_poudre
                T_voisin=[T(i+1,j) T(i,j+1) T(i+1,j+1) T(i-1,j) T(i,j-1) T(i-1,j-1) T(i-1,j+1) T(i+1,j-1)];
                x_voisin=[i+1 i i+1 i-1 i i-1 i-1 i+1];
                y_voisin=[j j+1 j+1 j j-1 j-1 j+1 j-1];
                [M,I]=min(T_voisin);
                if capa(x_voisin(I),y_voisin(I))==capa_air
                    T(x_voisin(I),y_voisin(I))=T(i,j);
                    capa(x_voisin(I),y_voisin(I))=capa(i,j);
                    diff(x_voisin(I),y_voisin(I))=diff(i,j);
                    capa(i,j)=capa_air;
                    diff(i,j)=diff_air;
                    Position_particules(i,j)=0;
                    Position_particules(x_voisin(I),y_voisin(I))=1;
                end
            end
        end
    end
    % Force de recul
    for i=2:distancex-1
        for j=2:distancey-1
            if capa(i,j)~=capa_air && T(i,j)<=10*T_fusion_poudre && T(i,j)>=T_fusion_poudre
                T_voisin=[T(i+1,j) T(i,j+1) T(i+1,j+1) T(i-1,j) T(i,j-1) T(i-1,j-1) T(i-1,j+1) T(i+1,j-1)];
                x_voisin=[1 0 1 -1 0 -1 -1 1];
                y_voisin=[0 1 1 0 -1 -1 1 -1];
                [M,I]=min(T_voisin);
                if capa(i-x_voisin(I),j-y_voisin(I))==capa_air
                    T(i-x_voisin(I),j-y_voisin(I))=T(i,j);
                    capa(i-x_voisin(I),j-y_voisin(I))=capa(i,j);
                    diff(i-x_voisin(I),j-y_voisin(I))=diff(i,j);
                    capa(i,j)=capa_air;
                    diff(i,j)=diff_air;
                    Position_particules(i,j)=0;
                    Position_particules(i-x_voisin(I),j-y_voisin(I))=1;
                end
            end
        end
    end
end