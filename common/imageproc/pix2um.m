function um = pix2um(params, pix)

%%
OX = params.ox;
OY = params.oy;
OZ = params.oz;
OO = [OX OY OZ]/1e3; % in um

SX = params.sx;
SY = params.sy;
SZ = params.sz;
OS = [SX SY SZ]/2^(params.level)/1e3;
%%
pix = pix-1+.5;
%%
um = pix.*(ones(size(pix,1),1)*OS) + ones(size(pix,1),1)*OO ;
