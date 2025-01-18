function Z=ZnmXY(X,Y,n,m)
% ZnmXY - évaluation d'une fonctions de base de Zernike EN CARTÉSIEN (normalisation écart-type à 1 sur disque rayon 1)
%       | computes a function of the Zernike basis in the CARTESIAN coordinate system (standard deviation normalisation to 1 on a disk of unit radius)
% USAGE: Z=ZnmXY(X,Y,n,m)
% ---
% X : tableau des X normalisés.
%   | array of normalized X.
% Y : tableau des Y normalisés. X et Y doivent avoir même taille (ou l'un des deux peut être scalaire)
%   | array of normalized Y. X and Y must have the same size (or either one can be scalar)
% n : scalaire entier positif ou nul indiquant le degré du polynôme de Zernike
%   | integer scalar specifying the degree of the Zernike polynomial
% m : scalaire entier indiquant l'indice 'angulaire' m du polynôme de Zernike
%     /!\ Il faut impérativement que 0 <= |m| <= n  et  m+n  pair! (Erreur à l'exécution sinon).
%   | integer scalar specifying the angular index of the Zernike polynomial
%   | /!\ m and n must be chosen such that 0 <= |m| <= n  and  m+n  even! (Error otherwise).       
% ---       
% Z   : tableau de même taille que X et Y donnant la valeur de  Znm(X,Y) = Rnm(rho).cos(m.theta).coef_normalisation
%     : pour m >= 0  et  Rnm(rho).sin(-m.theta).coef_normalisation pour m<0.
%     : Le coefficient de normalisation est calculé pour que l'écart-type sur le disque unité vaille 1 (et Z00 est constant à 1)
%     | Array of the same size as X and Y giving the value of  Znm(X,Y)=Rnm(rho).cos(m.theta).coef_normalization for m >= 0  and 
%     | Rnm(rho).sin(-m.theta).coef_normalisation for m<0.
%     | The normalization coefficient is computed such that the standard deviation on the unit circle equals 1 (and Z00 is constant to 1)     
% NOTA: Le lien entre X et Y et rho et theta est conforme aux usages de conception optique, c'est à dire
%     : # theta=0 dans la direction +y #. (On a donc y=rho.cos(theta) et x=-rho.sin(theta) ).
%     : Pour les systèmes centrés (i.e. à symétrie de révolution) avec champ suivant +y, les aberrations ne peuvent faire
%     : apparaître que des Znm avec m >=0.
% ¤¤¤¤| The relation between X and Y, and rho and theta is consistent with the convention used in optical design, that is to say
%     | # theta=0 in the direction +y #. (Thus  y=rho.cos(theta) and x=-rho.sin(theta) ).
%     | For a centered system (i.e. with a rotational symmetry) with the field along +y, aberrations can only lead to Znm polynomials
%     | with m >=0.
%
% Voir [Born&Wolf, Principles of Optics, (7th Ed, Cambridge 1999, ou Ed. précédente)] chapitre 9 pour la définition des polynômes 
% de Zernike... Avec les polynôme radiaux Rnm définis dans le [B&W], on a coef_normalisation=[2(n+1)/(1+kronecker(m,0))]^½ .
%| See [Born&Wolf, Principles of Optics, (7th Ed, Cambridge 1999, or previous Edition)] chapter 9 for the definition of the
%| Zernike polynomials... With the Rnm radial polynomials as defined in the [B&W], coef_normalization=[2(n+1)/(1+kronecker(m,0))]^½ .
% 
%  Utilise la fonction "foncZnm.m" de la BàO SupOptique
%| Uses the function "foncZnm.m" from the IOGS|‘SupOptique’ toolbox.
 
% © Institut d'Optique Graduate School - Hervé Sauer - 02 octobre 2012
% + H.S. - 09 décembre 2016 - Comments and Error messages in English

if ~isequal(size(X),size(Y)) && ~(isscalar(X) || isscalar(Y))
  error(['Les dimensions de X et Y doivent être identiques (ou l''un des deux doit être scalaire)' char(10) ...
          'The X and Y sizes should be identical (or one of X or Y argument must be scalar)'])
end

RHO=sqrt(X.^2+Y.^2);
TH=atan2(-X,Y); % place TH=0° suivant +y

Z=foncZnm(RHO,TH,n,m);
