// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
//
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
// Copyright IBM Corporation 1998, 2013. All Rights Reserved.
//
// Note to U.S. Government Users Restricted Rights:
// Use, duplication or disclosure restricted by GSA ADP Schedule
// Contract with IBM Corp.
// --------------------------------------------------------------------------

/****************************************************************** 
 * Life Game from Robert Bosch and Michael Trick, CP 2001, 
 * CPAIOR 2002.
 * Original model IP2 available at http://mat.gsia.cmu.edu/Life/
 * Basic integer program with birth constraints
 *
 * This model is greater than the size allowed in trial mode. 
 * You therefore need a commercial edition of CPLEX Studio to run this example. 
 * If you are a student or teacher, you can also get a full version through
 * the IBM Academic Initiative.
 *
******************************************************************/

main {
  var m1Source = new IloOplModelSource("lifegameip.mod");
  var def = new IloOplModelDefinition(m1Source);
  var cplex1 = new IloCplex();
  cplex1.readVMConfig(thisOplModel.resolvePath("process.vmc"));

  var opl1 = new IloOplModel(def, cplex1);
  opl1.generate();
  cplex1.solve();
  writeln(opl1.printSolution());

  if (cplex1.hasVMConfig()) {
    writeln("cplex1 has a VM config");
  } else {
    fail();
  }

  cplex1.delVMConfig();
  if (cplex1.hasVMConfig()) {
    writeln("cplex1 has a VM config");
  } else {
    writeln("cplex1 does not have anymore a VM config");
  }
}

