library(CNORode)
library(MEIGOR)
library(CellNOptR)

model = readSIF("rede1.sif")
cno_data=readMIDAS("exer2.csv")
cnolist=makeCNOlist(cno_data,subfield=FALSE)

plotModel(model,CNOlist =cnolist)




plotCNOlist(cnolist)

ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0.1,
                                   LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3,
                                   default_k = 0.5, default_tau = 1, opt_n = TRUE, opt_k = TRUE,
                                   opt_tau = TRUE, random = FALSE)

modelSim=plotLBodeModelSim(cnolist, model, ode_parameters,timeSignals=seq(0,2,0.5));

simulatedData=plotLBodeFitness(cnolist, model,
                               ode_parameters,
                               transfer_function = 3)

initial_pars=createLBodeContPars(model, LB_n = 1, LB_k = 0.1,
                                 LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, random = TRUE)

simulatedData=plotLBodeFitness(cnolist, model,initial_pars)
paramsGA = defaultParametersGA()
paramsGA$maxStepSize = 1
paramsGA$popSize = 50
paramsGA$iter = 100
paramsGA$transfer_function = 2
opt_pars=parEstimationLBode(cnolist,model,ode_parameters=initial_pars,
                            paramsGA=paramsGA)

simulatedData=plotLBodeFitness(cnolist, model,ode_parameters=opt_pars)                           
                           
requireNamespace("MEIGOR")
initial_pars=createLBodeContPars(model,
                                 LB_n = 1, LB_k = 0.1, LB_tau = 0.01, UB_n = 5,
                                 UB_k = 0.9, UB_tau = 10, random = TRUE)

fit_result_ess =
  parEstimationLBodeSSm(cnolist = cnolist,
                        model = model,
                        ode_parameters = initial_pars,
                        maxeval = 1e5,
                        maxtime = 20,
                        local_solver = "DHC",
                        transfer_function = 3
  )

                           
simulatedData=plotLBodeFitness(cnolist, model,
                               ode_parameters=fit_result_ess,
                               transfer_function = 3)                          



