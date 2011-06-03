optParamatersDEoptim <-
function(pars,ub,lb,f)
{  
   library(DEoptim);
   DEoptim(fn =f, lower =lb , upper =ub, control = list(NP = 10, itermax = 10, trace = TRUE));
}

