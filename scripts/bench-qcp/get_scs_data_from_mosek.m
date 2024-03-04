function [scs_data, scs_cones] = get_scs_data_from_mosek(mosek_res)

    res = mosek_res;
    
    if res.rcodestr ~= 'MSK_RES_OK'
        error("mosek read data error");
    end
    
     if res.prob.cones.type(1) == 1
        error("RSOC not supported by scs\n");
    end
    
    soc = [];
    
    if(length(res.prob.cones.subptr) > 1)
        for ii = 1:length(res.prob.cones.subptr) - 1
            soc = [soc, res.prob.cones.subptr(ii + 1) - res.prob.cones.subptr(ii)];
        end
    end
    
    soc = [soc, length(res.prob.cones.sub) - res.prob.cones.subptr(end) + 1];
        
    
    ind_blc = ~isinf(res.prob.blc);
    ind_buc = ~isinf(res.prob.buc);
    ind_blx = ~isinf(res.prob.blx);
    ind_bux = ~isinf(res.prob.bux);
    
    blc = res.prob.blc(ind_blc);
    buc = res.prob.buc(ind_buc);
    blx = res.prob.blx(ind_blx);
    bux = res.prob.bux(ind_bux);
    
    len_blc = length(blc);
    len_buc = length(buc);
    len_blx = length(blx);
    len_bux = length(bux);
    
    A = res.prob.a;
    c = res.prob.c;
    
    [m,n] = size(A);
    %
    A_new = [ [A(ind_blc,:);sparse(1:len_blx,find(ind_blx==1),ones(len_blx,1),len_blx,n)], -speye(len_blc + len_blx), sparse(len_blc + len_blx,len_buc + len_bux);
              [A(ind_buc,:);sparse(1:len_bux,find(ind_bux==1),ones(len_bux,1),len_bux,n)], sparse(len_buc + len_bux,len_blc + len_blx), speye(len_buc + len_bux)];
          
    b_new = [blc;blx;buc;bux];
    
    c_new = [c;zeros(len_blc + len_blx + len_buc + len_bux, 1)];
    
    free_cone = setdiff([1:n], res.prob.cones.sub);
    
    A_new = [A_new(:,free_cone), A_new(:,n+1:end), A_new(:,res.prob.cones.sub)];
    c_new = [c_new(free_cone); c_new(n+1:end); c_new(res.prob.cones.sub)];
    
    scs_data.A = A_new';
    scs_data.c = -b_new;
    
    if res.prob.objsense == 'minimize'
        scs_data.b = c_new;
    else
        scs_data.b = -c_new;
    end
        
    scs_cones.z = length(free_cone);
    scs_cones.l = len_blc + len_blx + len_buc + len_bux;
    scs_cones.q = soc;
    
end