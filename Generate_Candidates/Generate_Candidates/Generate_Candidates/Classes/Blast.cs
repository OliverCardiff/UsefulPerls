using System;

namespace Generate_Candidates
{
	public class Blast
	{
		public Scaffold _A;
		public Scaffold _B;

		public int _st;
		public int _ed;
		public int _tSt;
		public int _tEd;
		public double _percID;

		public Blast (Scaffold A, Scaffold B, double percID, int st, int ed, int tst, int ted)
		{
			_percID = percID;
			_A = A;
			_B = B;
			_st = st;
			_ed = ed;
			_tSt = tst;
			_tEd = ted;
		}

		public bool Overlaps(Allele oth)
		{
			if ((_st <= oth._ed && _st >= oth._st) || (_ed <= oth._ed && _ed >= oth._st)) {
				return true;
			} else {
				return false;
			}
		}
	}
}

