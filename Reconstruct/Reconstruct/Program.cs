using System;

namespace Reconstruct
{
	class MainClass
	{

		public static int LengthDiv = 5;

		public static void Main (string[] args)
		{
			string mode = "";
			string gen_file = "";
			string links_file = "";
			string out_file = "";

			try{
			mode = args [0];
			gen_file = args[1];
			links_file = args [2];
			out_file = args [3];
			} catch(Exception e) {
				Console.WriteLine ("Input Error!");
				Console.WriteLine (e.Message);
				Console.WriteLine ("Usage: Reconstruct <mode> <genome> <link_file> <output>");
				Console.WriteLine ("Mode Options: 'scaffold' 'overlap'");
				return;
			}
				
			Console.WriteLine ("Reading .Fai index..");
			Genome genome = new Genome (gen_file, links_file, out_file);
			genome.LiteLoadGenome ();
			Console.WriteLine ("Reading Links File...");
			genome.LoadAlleles ();
			if (mode.Equals ("overlap")) {
				Console.WriteLine ("Ran in 'overlap' mode..");
				Console.WriteLine ("Eliminating small fragments..\n");
				//genome.ExpandReflections (false);
				//genome.MergeOverlaps ();
				//genome.ExpandReflections (false);
				//genome.MergeOverlaps ();
				//genome.EliminateSmallFrags ();
				//genome.ExpandReflections (false);
				//genome.MergeOverlaps ();
				genome.EliminateSmallFrags ();

				genome.KeepJoins ();
				genome.BuildMetaScaffolds ();

				Console.WriteLine ("Reading in Genome: ..");
				genome.LoadGenome ();

				Console.WriteLine ("Writing out everything...");
				genome.WriteGenome (true);

			} else if (mode.Equals ("scaffold")) {
				Console.WriteLine ("Ran in 'scaffold' mode..");
				Console.WriteLine ("Making hypercontigs..");
				//genome.build!
			} else {
				Console.WriteLine ("ERROR: Enter 'scaffold' or 'overlap' as the first argument\n");
			}
		}
	}
}
