using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Reconstruct
{
	public class AlleleGroup
	{
		List<ScaffIndices> _indicies;
		List<ScaffIndices> _coords;
		List<Scaffold> _scaffs;
		List<int> _Npads;

		public AlleleGroup()
		{
			_indicies = new List<ScaffIndices> ();
			_scaffs = new List<Scaffold> ();
			_Npads = new List<int> ();
			_coords = new List<ScaffIndices> ();
		}
	}
	public class HyperContig
	{
		public Scaffold _starter;

		protected AlleleGroup _A;
		protected AlleleGroup _B;

		public HyperContig (Scaffold st)
		{
			_starter = st;
			_A = new AlleleGroup ();
			_B = new AlleleGroup ();
		}

		public int BuildPairs()
		{
			return 0;
		}

	}
}

