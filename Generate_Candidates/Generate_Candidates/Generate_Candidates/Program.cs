using System;

//Arguments Genome Wig Blast Prefix

namespace Generate_Candidates
{
	class MainClass
	{
		public static int LOWCOV = 60;
		public static int HIGHCOV = 110;
		public static int RANGE = 15;
		public static int BLASTPERCID = 99;

		public static void Main (string[] args)
		{
			string gen_file = args[0];
			string wig_file = args [1];
			string blast_file = args [2];
			string prefix = args [3];

			Console.WriteLine ("Reading Genome..");
			Genome genome = new Genome (gen_file, wig_file, blast_file);
			genome.LoadGenome ();
			Console.WriteLine ("Reading Blast file..");
			genome.LoadBlast ();
			Console.WriteLine ("Reading Wig..");
			genome.LoadWigs ();
			Console.WriteLine ("Generating all non-collapse alleles\n");
			genome.GenerateAllAlleles ();
			Console.WriteLine ("Piling up Blast Results on Alleles\n");
			genome.GenerateAllAnchors ();
			Console.WriteLine ("Writing out anchor file to: " + prefix + "..\n");
			genome.WriteLinksFile (prefix);
			Console.WriteLine ("Done");
		}
	}
}
