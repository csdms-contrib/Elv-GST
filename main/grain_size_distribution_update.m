%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 81 $
%$Date: 2017-04-20 15:51:02 +0200 (Thu, 20 Apr 2017) $
%$Author: V $
%$Id: grain_size_distribution_update.m 81 2017-04-20 13:51:02Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/grain_size_distribution_update.m $
%
%grain_size_distribution_update updates the mass at the active layer and the substrate and the substrate thickness
%
%[Mak_new,msk_new,Ls_new,La_ne,etab_new,ell_idx,out_en,pmm]=grain_size_distribution_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,bc,u,h,Cf,input,fid_log,kt)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.
%
%160429
%   -V. Bug solved. In case the active layer thickness is updated the
%   derivatives need to be computed with the old one because is the one
%   used in the computation of the sediment transport (nt changed in
%   Hirano)
%
%166614
%   -V. Output pmm for adjusting CFL condition
%
%160701
%   -V. Hoey and Ferguson

function [Mak_new,msk_new,Ls_new,La_ne,etab_new,ell_idx,out_en,pmm]=grain_size_distribution_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,bc,u,h,Cf,input,fid_log,kt)
%comment out fot improved performance if the version is clear from github
% version='4';
% if kt==1; fprintf(fid_log,'grain_size_distribution_update version: %s\n',version); end 

%%
%% RENAME
%%

%%
%% GSD UPDATE
%%

switch input.mor.gsdupdate
%% NO UPDATE    
    case 0
        Mak_new=Mak;
        msk_new=msk;
        Ls_new=Ls;
        ell_idx=false(1,input.mdv.nx);
        La_ne=La; %non-elliptic La
        out_en=NaN;
        etab_new=etab;
        pmm.alpha=ones(1,input.mdv.nx);
        pmm.beta=ones(1,input.mdv.nx);
        
%% HIRANO        
    case 1       
        [msk_new,Ls_new,fIk,detaLa]=substrate_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,input,fid_log,kt);
        Mak_new=active_layer_mass_update(Mak,detaLa,fIk,qbk,bc,input,fid_log,kt);
        if input.mor.ellcheck
            [ell_idx,out_en]=elliptic_nodes(u,h,Cf,La_old,qbk,Mak,fIk,input,fid_log,kt);
        else
            ell_idx=false(1,input.mdv.nx);
            out_en=NaN;
        end
        La_ne=La; %non-elliptic La
        etab_new=etab;
        pmm.alpha=ones(1,input.mdv.nx);
        pmm.beta=ones(1,input.mdv.nx);
        
%% ELI 1
    case {2,3}
        [~,~,fIk,~]=substrate_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,input,fid_log,kt); %substrate to obtain fIk
        [ell_idx,out_en]=elliptic_nodes(u,h,Cf,La_old,qbk,Mak,fIk,input,fid_log,kt); %La_old because the derivatives are of the transport before updating La
        La_ne=eli_1(La,ell_idx,out_en,input,fid_log,kt);
        [msk_new,Ls_new,fIk,detaLa]=substrate_update(Mak,msk,Ls,La_old,La_ne,etab_old,etab,qbk,input,fid_log,kt); %obtain the new substrate
        Mak_new=active_layer_mass_update(Mak,detaLa,fIk,qbk,bc,input,fid_log,kt); %obtain the new mass at the active layer
        etab_new=etab;
        pmm.alpha=ones(1,input.mdv.nx);
        pmm.beta=ones(1,input.mdv.nx);
        
%% PRE-CONDITIONING MASS MATRIX
    case {4,5,6,7,8}
        [~,~,fIk,~]=substrate_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,input,fid_log,kt); %substrate to obtain fIk
        [ell_idx,out_en]=elliptic_nodes(u,h,Cf,La_old,qbk,Mak,fIk,input,fid_log,kt); %La_old because the derivatives are of the transport before updating La
        pmm=preconditioning_mass_matrix(ell_idx,out_en,u,input,fid_log,kt);
        etab_new=bed_level_update_pmm(etab_old,qbk,bc,pmm,input,fid_log,kt); %corrected etab
        [msk_new,Ls_new,fIk,detaLa]=substrate_update(Mak,msk,Ls,La_old,La,etab_old,etab_new,qbk,input,fid_log,kt); %obtain the new substrate
        Mak_new=active_layer_mass_update_pmm(Mak,detaLa,fIk,qbk,bc,pmm,input,fid_log,kt); %obtain the new mass at the active layer        
        La_ne=La; %non-elliptic La
end %gsdupdate


end %function
    
    
    




