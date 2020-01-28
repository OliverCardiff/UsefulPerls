using System;
using System.Collections.Generic;
using System.IO;

namespace Generate_Candidates
{
	public class Scaffold
	{
		/*public static bool operator == (Scaffold a, Scaffold b)
		{
			if (a._name.Equals (b._name)) {
				return true;
			} else {
				return false;
			}
		}
		public static bool operator != (Scaffold a, Scaffold b)
		{
			if (a._name.Equals (b._name)) {
				return false;
			} else {
				return true;
			}
		}
		public override bool Equals (object obj)
		{
			if (obj.ToString().Equals (_name)) {
				return true;
			} else {
				return false;
			}
		}*/
		public override string ToString ()
		{
			return _name;
		}
		protected List<int> _wig;
		protected List<Blast> _blasts;
		protected List<Allele> _alleles;

		public string _name;
		public int _length;
		public string _sequence;

		public Scaffold(string id, ref string seq)
		{
			_sequence = seq;
			_name = id;
			_length = seq.Length;
			_blasts = new List<Blast> ();
			_alleles = new List<Allele> ();
		}

		public void GiveWig(List<int> wig)
		{
			_wig = wig;
		}

		public void GiveBlast(Blast bla)
		{
			_blasts.Add(bla);
		}

		public int PileupBlasts()
		{
			int cnt = 0;
			foreach (Blast b in _blasts) {
				foreach (Allele a in _alleles) {

					if (b.Overlaps (a)) {
						a.AddBlast (b);
					}
				}
			}
			foreach (Allele a in _alleles) {
				cnt += a.GenerateAnchorBases ();
			}
			return cnt;
		}

		public void WriteLinks(StreamWriter sw)
		{
			foreach (Allele a in _alleles) {
				a.WriteSpans (sw);
			}
		}

		public void GenerateAlleles()
		{
			int primaryBroad = 0;
			int primaryNarrow = 0;
			int finalMerge = 0;
			int joins = 0;
			int total = 0;
			int tot_len = 0;
			List<int> broadWindow = new List<int>();
			List<int> narrowWindow = new List<int>();

			RunWindows (broadWindow, narrowWindow);

			List<Allele> nAllele = PopulateAlleleList(narrowWindow, 100);
			List<Allele> bAllele = PopulateAlleleList(broadWindow, 500);

			primaryBroad = bAllele.Count;
			primaryNarrow = nAllele.Count;

			List<Allele> Merged = MergeAlleleLists (bAllele, nAllele);

			finalMerge = Merged.Count;

			_alleles = JoinAdjacentAlleles (Merged);

			total = _alleles.Count;
			joins = finalMerge - total;
			tot_len = SumofAlleles (_alleles);

			if (total > 99999) {
				Console.WriteLine ("Allele Generation Report for: " + _name);
				Console.WriteLine ("Inital Broad:\t" + primaryBroad);
				Console.WriteLine ("Inital Narrow:\t" + primaryNarrow);
				Console.WriteLine ("Merged Count:\t" + finalMerge);
				Console.WriteLine ("Joins Made:\t" + joins);
				Console.WriteLine ("Final Count:\t" + total);
				int cnt2 = 0;
				foreach (Allele a in _alleles) {
					cnt2++;
					if (a.Length > 2000) {
						Console.WriteLine ("Allele " + cnt2 + " st: " + a._st + "\ted: " + a._ed + "\tlen: " + a.Length);
					}
				}
				Console.WriteLine ("Final Sum:\t" + tot_len);
				Console.ReadLine ();
			}
		}

		protected int SumofAlleles(List<Allele> al)
		{
			int accu = 0;
			foreach (Allele a in al) {
				accu += a.Length;
			}
			return accu;
		}

		protected List<Allele> JoinAdjacentAlleles(List<Allele> all)
		{
			int len = all.Count;
			List<Allele> final = new List<Allele> ();

			for (int i = 0; i < len - 1; i++) {
				for (int j = i + 1; j < len; j++) {
					int prox = all [i].GetProximity (all [j]);
					int lenBoth = Math.Max(all[i].Length, all[j].Length);
					int bound = Math.Min (lenBoth / 10, 2000);

					if (prox < bound) {
						all [j] = new Allele (all [i], all [j]);
						all [i] = null;
						break;
					}
				}
			}

			foreach (Allele a in all) {
				if (a != null) {
					final.Add (a);
				}
			}

			return final;
		}

		protected List<Allele> MergeAlleleLists(List<Allele> broad, List<Allele> narrow)
		{
			int nLen = narrow.Count;
			int bLen = broad.Count;

			for (int j = 0; j < bLen; j++) 
			{
				for (int i = 0; i < nLen; i++) {
					if (narrow [i] != null) {
						if (broad [j].Contains (narrow [i])) {
							narrow [i] = null;
						} else if (broad [j].Overlaps (narrow [i])) {
							broad [j] = new Allele (broad [j], narrow [i]);
							narrow [i] = null;
						}
					}
				}
			}

			return broad;
		}

		protected List<Allele> PopulateAlleleList(List<int> wiggo, int minAl)
		{
			List<Allele> nAllele = new List<Allele> ();

			bool current = false;
			int tolerance = 3;
			bool inTolerance = false;
			int st = 0; int ed = 0;

			for (int i = 0; i < wiggo.Count; i++) {

				if (wiggo [i] < MainClass.LOWCOV + MainClass.RANGE) {
					if (!current) {
						st = i * 10;
						current = true;
					}
					inTolerance = false;
					tolerance = Math.Min (tolerance + 1, 3);

				} else {
					
					if (current) {
						if (tolerance > 0) {
							if (inTolerance == false) {
								ed = i * 10;
							}
							inTolerance = true;
							tolerance--;
						} else {
							if (ed - st > minAl) {
								nAllele.Add (new Allele (this, st, ed));
							}
							current = false;
							st = 0; ed = 0;
							tolerance = 3;
							inTolerance = false;
						}
					}
				}
			}
			return nAllele;
		}

		protected void RunWindows(List<int> broadWindow, List<int> narrowWindow)
		{
			if (_wig == null) {
				FillNullWig ();
			}
			int len = _wig.Count;

			for (int i = 0; i < len; i++) {

				double cnt = 0; double accu = 0;

				int winLenAB = Math.Max (i - 15, 0);
				int winLenBB = Math.Min (i + 15, len - 1);

				int winLenAN = Math.Max (i - 5, 0);
				int winLenBN = Math.Min (i + 5, len - 1);

				InnerLoop (i, winLenAB, ref accu, ref cnt);
				InnerLoop (i, winLenBB, ref accu, ref cnt);

				broadWindow.Add((int)(accu / cnt));

				accu = 0; cnt = 0;

				InnerLoop (i, winLenAN, ref accu, ref cnt);
				InnerLoop (i, winLenBN, ref accu, ref cnt);

				narrowWindow.Add((int)(accu / cnt));
			}
		}

		protected void FillNullWig()
		{
			int len = _length / 10;

			_wig = new List<int> ();
			for (int i = 0; i < len; i++) {
				_wig.Add (MainClass.HIGHCOV);
			}
		}

		protected void InnerLoop(int i, int winLen, ref double accu, ref double cnt)
		{
			int increment = -1;
			if (i < winLen) {
				increment = 1;
			}
			winLen += increment;

			for (int j = i; j != winLen; j += increment) {
				cnt++;
				Accrue(MainClass.HIGHCOV, ref accu, _wig[j]);
			}
		}

		protected void Accrue(int bound, ref double accu, int lev)
		{
			if (lev == 0) {
				accu += bound;
			} else {
				accu += Math.Min (lev, bound);
			}
		}
	}
}

