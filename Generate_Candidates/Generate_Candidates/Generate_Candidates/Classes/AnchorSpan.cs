using System;
using System.Collections.Generic;
using System.IO;

namespace Generate_Candidates
{
	public struct SyntenyAnchor
	{
		public SyntenyAnchor(SyntenyAnchor oth)
		{
			Local = oth.Local;
			Target = oth.Target;
			Score = oth.Score;
		}
		public SyntenyAnchor(int st, int ed, int score)
		{
			Local = st;
			Target = ed;
			Score = score;
		}
		public int Local;
		public int Target;
		public int Score;
	};
	public class AnchorSpan
	{
		public int _st;
		public int _ed;

		public int _tSt;
		public int _tEd;
		List<AnchorBase> _bases;
		public bool _reversed;
		public SyntenyAnchor _anchor;
		public Scaffold _target;

		public int Length { get { return _ed - _st; } }
		public int Count { get { return _bases.Count; } }

		public AnchorSpan (List<AnchorBase> bases)
		{
			_target = bases [0]._target;
			_bases = bases;
			SetLimits ();
			SetReversed ();
			_anchor = GetSyntenyAnchorPoint ();
		}

		public void Write(StreamWriter sw, Scaffold parent)
		{
			sw.WriteLine (parent._name + "\t" + _target._name + "\t" + _st + "\t" + _ed + "\t" + _tSt + "\t" + _tEd);
		}

		public bool Overlaps(AnchorSpan oth)
		{
			if ((_st <= oth._ed && _st >= oth._st) || (_ed <= oth._ed && _ed >= oth._st)) {
				return true;
			} else {
				return false;
			}
		}

		public bool Contains(AnchorSpan oth)
		{
			if (_st >= oth._st && _ed <= oth._ed) {
				return true;
			} else {
				return false;
			}
		}

		public void DestroyOverlap(int st, int ed)
		{
			List<AnchorBase> nextList = new List<AnchorBase> ();

			foreach (AnchorBase ab in _bases) {
				if (!ab.IsWithin (st, ed)) {
					nextList.Add (ab);
				}
			}

			_bases = nextList;
			SetLimits ();
			_anchor = GetSyntenyAnchorPoint ();
		}

		protected void SetLimits()
		{
			_st = _bases [0]._loc;
			_ed = _bases [_bases.Count - 1]._loc;

			int min = 999999999;
			int max = 0;

			foreach (AnchorBase ab in _bases) {
				if (ab._tarLoc < min) {
					min = ab._tarLoc;
				}
				if (ab._tarLoc > max) {
					max = ab._tarLoc;
				}
			}
			if (_reversed) {
				_tSt = max;
				_tEd = min;
			} else {
				_tSt = min;
				_tEd = max;
			}
		}

		protected SyntenyAnchor GetSyntenyAnchorPoint()
		{
			int quart = _bases.Count / 4;
			int median = _bases.Count / 2;
			int quart3 = median + quart;

			int medSt = _bases [median]._loc;
			int medEd = _bases [median]._tarLoc;

			int qrtSt = _bases [quart]._loc;
			int qrtEd = _bases [quart]._tarLoc;

			int q3St = _bases [quart3]._loc;
			int q3Ed = _bases [quart3]._tarLoc;

			SyntenyAnchor Q1 = OptimiseAnchor (new SyntenyAnchor (qrtSt, qrtEd, 0), 50, 5);
			SyntenyAnchor Med = OptimiseAnchor (new SyntenyAnchor (medSt, medEd, 0), 50, 5);
			SyntenyAnchor Q3 = OptimiseAnchor (new SyntenyAnchor (q3St, q3Ed, 0), 50, 5);

			SyntenyAnchor chosen;

			if (Q1.Score < Med.Score && Q1.Score < Q3.Score) {
				chosen = Q1;
			} else if (Med.Score < Q3.Score) {
				chosen = Med;
			} else {
				chosen = Q3;
			}

			SyntenyAnchor Opt1 = OptimiseAnchor (chosen, 10, 5);
			SyntenyAnchor Opt2 = OptimiseAnchor (Opt1, 1, 10);

			return Opt2;
		}

		protected SyntenyAnchor OptimiseAnchor(SyntenyAnchor syn, int increment, int steps)
		{
			int st = syn.Local;
			int ed = syn.Target;

			ed -= increment * steps;

			int lowestSum = 999999999;
			int edVal = 0;

			for (int i = 0; i <= steps * 2; i++) {

				int test = Math.Abs(SyntenySum (st, ed));
				if (test < lowestSum) {
					lowestSum = test;
					edVal = ed;
				}
				ed += increment;
			}

			return new SyntenyAnchor (st, edVal, lowestSum);
		}

		protected void SetReversed()
		{
			int revCnt = 0;
			int normCnt = 0;
			foreach (AnchorBase ab in _bases) {
				if (ab._reversed) {
					revCnt++;
				} else {
					normCnt++;
				}
			}
			if (revCnt > normCnt) {
				_reversed = true;
			} else {
				_reversed = false;
			}
		}

		protected int SyntenySum(int stA, int stB)
		{
			int diff = 0;

			foreach (AnchorBase ab in _bases) {
				if (_reversed) {
					diff += (ab._loc - stA) - (stB - ab._tarLoc);
				} else {
					diff += (ab._loc - stA) - (ab._tarLoc - stB);
				}
			}
			return diff;
		}

		public double GetDensitySpan(int st, int ed)
		{
			double len = (double)(ed - st);
			double accu = 0;

			foreach (AnchorBase ab in _bases) {
				if (ab.IsWithin (st, ed)) {
					accu++;
				}
			}

			return accu / len;
		}
	}
}

