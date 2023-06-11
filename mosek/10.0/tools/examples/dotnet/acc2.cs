/*
  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File:      acc2.cs

  Purpose :   Tutorial example for affine conic constraints.
              Models the problem:
 
              maximize c^T x
              subject to  sum(x) = 1
                          gamma >= |Gx+h|_2

              This version inputs the linear constraint as an affine conic constraint.
*/
using System;

namespace mosek.example
{
  class msgclass : mosek.Stream
  {
    string prefix;
    public msgclass (string prfx)
    {
      prefix = prfx;
    }

    public override void streamCB (string msg)
    {
      Console.Write ("{0}{1}", prefix, msg);
    }
  }

  public class acc1
  {
    public static void Main ()
    {
      /* Problem dimensions */
      const int n = 3;
      const int k = 2;

      int i,j;
      long quadDom, zeroDom;

      // Since the value infinity is never used, we define
      // 'infinity' symbolic purposes only
      double infinity = 0;

      // Create a task object.
      using (mosek.Task task = new mosek.Task()) {
        // Directs the log task stream to the user specified
        // method msgclass.streamCB
        task.set_Stream (mosek.streamtype.log, new msgclass (""));

        // Create n free variables
        task.appendvars(n);
        task.putvarboundsliceconst(0, n, mosek.boundkey.fr, -infinity, infinity);

        // Set up the objective
        double[] c = {2, 3, -1};
        int[] cind = {0, 1, 2};  
        task.putobjsense(mosek.objsense.maximize);
        task.putclist(cind, c);

        // Set AFE rows representing the linear constraint
        task.appendafes(1);
        task.putafeg(0, -1.0);
        for(i = 0; i < n; i++) task.putafefentry(0, i, 1.0);

        // F matix in sparse form
        long[]   Fsubi = {2, 2, 3, 3};   // The G matrix starts in F from row 2
        int[]    Fsubj = {0, 1, 0, 2};
        double[] Fval  = {1.5, 0.1, 0.3, 2.1};
        // Other data
        double[] h     = {0, 0.1};
        double   gamma = 0.03;

        task.appendafes(k + 1);
        task.putafefentrylist(Fsubi, Fsubj, Fval);
        task.putafeg(1, gamma);
        task.putafegslice(2, k+2, h);

        // Define domains
        zeroDom = task.appendrzerodomain(1);
        quadDom = task.appendquadraticconedomain(k + 1);

        // Create the linear ACC
        long[] afeidxZero = {0};
        task.appendacc(zeroDom,    // Domain index
                       afeidxZero,// Indices of AFE rows
                       null);      // Ignored

        // Create the quadratic ACC
        long[] afeidxQuad = {1, 2, 3};
        task.appendacc(quadDom,    // Domain index
                       afeidxQuad, // Indices of AFE rows
                       null);      // Ignored

        // Solve the problem
        task.optimize();

        // Print a summary containing information
        // about the solution for debugging purposes
        task.solutionsummary(mosek.streamtype.msg);

        /* Get status information about the solution */
        mosek.solsta solsta;
        task.getsolsta(mosek.soltype.itr, out solsta);

        switch (solsta)
        {
          case mosek.solsta.optimal:
            // Fetch solution
            double[] xx  = new double[n];
            task.getxx(mosek.soltype.itr, // Interior solution.
                       xx);
            Console.WriteLine ("Optimal primal solution");
            for (j = 0; j < n; ++j)
              Console.WriteLine ("x[{0}]: {1}", j, xx[j]);

            // Fetch doty dual of the ACC
            double[] doty  = new double[k+1];
            task.getaccdoty(mosek.soltype.itr, // Interior solution.
                            1,                 // ACC index
                            doty);
            Console.WriteLine ("Dual doty of ACC");
            for (j = 0; j < k+1; ++j)
              Console.WriteLine ("doty[{0}]: {1}", j, doty[j]);

            // Fetch activity of the ACC
            double[] activity  = new double[k+1];
            task.evaluateacc(mosek.soltype.itr, // Interior solution.
                             1,                 // ACC index
                             activity);
            Console.WriteLine ("Activity of ACC");
            for (j = 0; j < n; ++j)
              Console.WriteLine ("activity[{0}]: {1}", j, activity[j]);
            break;

          case mosek.solsta.dual_infeas_cer:
          case mosek.solsta.prim_infeas_cer:
            Console.WriteLine("Primal or dual infeasibility.\n");
            break;
          case mosek.solsta.unknown:
            Console.WriteLine("Unknown solution status.\n");
            break;
          default:
            Console.WriteLine("Other solution status");
            break;
        }
      }
    }
  }
}
