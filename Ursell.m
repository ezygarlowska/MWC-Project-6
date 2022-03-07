%this function computes the Ursell number from k, h and Hrms
function [Ur] = Ursell(k,h,Hrms)
    a=0.5*sqrt(2).*Hrms;
    Ur=(3*a.*k)./(4*(k.*h).^3);
end
