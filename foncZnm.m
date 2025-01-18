function fZ=foncZnm(u,th,n,m)
% foncZnm: évaluation des fonctions de base de Zernike en polaire (normalisation écart-type à 1 sur disque rayon 1)
% USAGE: fZ=foncZnm(u,th,n,m)
% u : tableau des rayons normalisés
% th: tableau de même taille que u , les angles [en radians]
% n : scalaire entier positif indiquant l'indice 'radiale' n du polynôme de Zernike
% m : scalaire entier indiquant l'indice 'angulaire' m du polynôme de Zernike
%     /!\ Il faut impérativement que 0 <= |m| <= n  et  m+n  paire! (Erreur à l'éxécution sinon).
%        
% fZ: tableau de même taille que u et th donnant la valeur de  Rnm(u).cos(m.theta).coef_normalisation
%     pour m >= 0  et  Rnm(u).sin(-m.theta).coef_normalisation pour m<0.
%
% Voir [Born&Wolf, Principles of Optics, (7th Ed, Cambridge 1999, ou Ed. précédente)] chapitre 9 pour la
% définition des polynômes de Zernike...
% coef_normalisation vaut [2(n+1)/(1+kronecker(m,0))]^½ 
% NOTA: Les polynômes de Zernike ne sont DÉFINIS que pour n positif ou 0, |m| <= n
% ^^^^  et  m  de même parité que  n .
% ====   m  positif (>=0) donne un terme angulaire en cos(m.th),
% ====   m  négatif (< 0) donne un terme angulaire en sin(-m.th)
%           {Les m<0 ne sont jamais présents pour les systèmes centrés avec champ en y (theta=0 suivant +y)}
 
% © Institut d'Optique Graduate School - TD 2A, Conception de Systèmes Optiques.
% Hervé Sauer - 29/10/2011 - D'après version antérieure TD Ab&Diff2A, cours Pierre Chavel.
% H.S. - 09/12/2016 - Erreurs et warnings bilingues Fr & Ang...


if ~isequal(size(u),size(th))
  error(['les arguments  u  et  th  doivent être des tableaux de mêmes dimensions!' char(10) ...
         'Arguments u and th must be array of the same sizes'])
end

if ~isscalar(n) || ~isscalar(m) || mod(n,1)~=0 || mod(m,1)~=0
  error(['les indices  n  et  m  doivent être des scalaires entiers!' char(10)...
         ' Arguments  n  and  m  must be integer scalars'])
end

if n<0 || abs(m)>n || mod(n+m,2)~=0
  error(['Il faut impérativement  0<= |m| <= n  et  n+m paire...' char(10) ...
         'It is mandatory that  0<= |m| <= n  and  n+m be even'] )
end

if n>20
  warning(['/!\ Le degré du polynôme de Zernike est TRÈÈÈS élevé! (Calcul sensé?)' char(10) ...
           '/!\ The Zernike polynomial degree is VERY high! (Is the computation really meaningful?)'])
end


% Calcul des coefficients du polynôme radiale Rnm(u):
% Voir [Born&Wolf, Principles of Optics, 6th Ed, Pergamon 1980], §9.2, p.465...
mz=abs(m);
CoefRnm=zeros(1,n+1);
s=0:(n-mz)/2;
tmp=fact([n-s;s;(n+mz)/2-s;(n-mz)/2-s]); % fact: BàO SupOptique (factoriel vectorisé)
CoefRnm(n+1-(n-2*s))=(-1).^s .* tmp(1,:)./prod(tmp(2:end,:));


% Évaluation de la fonction de base de Zernike en (u,th):
if m>=0
  fZ=polyval(CoefRnm,u).*cos(m*th);
else
  fZ=polyval(CoefRnm,u).*sin(-m*th);
end


fZ=fZ*sqrt(2*(n+1)/(1+double(m==0))); % NORMALISATION pour écart-type à 1 sur le disque de rayon unité

