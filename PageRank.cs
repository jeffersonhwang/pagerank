public class PageRank
{
	#region Private Fields

	ArrayList _incomingLinks, _leafNodes;
	Vector<double> _numLinks;
	double _alpha, _convergence;
	int _checkSteps;

	#endregion

	#region Constructor

	public PageRank(ArrayList linkMatrix, double alpha = 0.85, double convergence = 0.0001, int checkSteps = 10)
	{
		Tuple<ArrayList, Vector<double>, ArrayList> tuple = TransposeLinkMatrix(linkMatrix);
		_incomingLinks = tuple.Item1;
		_numLinks = tuple.Item2;
		_leafNodes = tuple.Item3;
		_alpha = alpha;
		_convergence = convergence;
		_checkSteps = checkSteps;
	}

	#endregion

	#region Methods

	/// <summary>
	/// Convenience wrap for the link matrix transpose and the generator.
	/// See PageRankGenerator method for parameter descriptions
	/// </summary>
	public double[] ComputePageRank()
	{
		Vector<double> final = null;
		foreach (Vector<double> generator in PageRankGenerator(_incomingLinks, _numLinks, _leafNodes, _alpha, _convergence, _checkSteps))
			final = generator;

		return final.ToArray();
	}

	/// <summary>
	/// Transposes the link matrix which contains the links from each page. 
	/// Returns a Tuple of:  
	/// 1) pages pointing to a given page, 
	/// 2) how many links each page contains, and
	/// 3) which pages contain no links at all. 
	/// We want to know is which pages
	/// </summary>
	/// <param name="outGoingLinks">outGoingLinks[i] contains the indices of the pages pointed to by page i</param>
	/// <returns>A tuple of (incomingLinks, numOutGoingLinks, leafNodes)</returns>
	protected Tuple<ArrayList, Vector<double>, ArrayList> TransposeLinkMatrix(ArrayList outGoingLinks)
	{
		int nPages = outGoingLinks.Count;

		// incomingLinks[i] will contain the indices jj of the pages
		// linking to page i
		ArrayList incomingLinks = new ArrayList(nPages);
		for (int i = 0; i < nPages; i++)
			incomingLinks.Add(new List<int>());

		// the number of links in each page
		Vector<double> numLinks = new DenseVector(nPages);

		// the indices of the leaf nodes
		ArrayList leafNodes = new ArrayList();
		for (int i = 0; i < nPages; i++)
		{
			List<int> values = outGoingLinks[i] as List<int>;
			if (values.Count == 0)
				leafNodes.Add(i);
			else
			{
				numLinks[i] = values.Count;
				// transpose the link matrix
				foreach (int j in values)
				{
					List<int> list = (List<int>)incomingLinks[j];
					list.Add(i);
					incomingLinks[j] = list;
				}
			}
		}

		return new Tuple<ArrayList, Vector<double>, ArrayList>(incomingLinks, numLinks, leafNodes);
	}

	/// <summary>
	/// Computes an approximate page rank vector of N pages to within some convergence factor.
	/// </summary>
	/// <param name="at">At a sparse square matrix with N rows. At[i] contains the indices of pages jj linking to i</param>
	/// <param name="leafNodes">contains the indices of pages without links</param>
	/// <param name="numLinks">iNumLinks[i] is the number of links going out from i.</param>
	/// <param name="alpha">a value between 0 and 1. Determines the relative importance of "stochastic" links.</param>
	/// <param name="convergence">a relative convergence criterion. Smaller means better, but more expensive.</param>
	/// <param name="checkSteps">check for convergence after so many steps</param>
	protected IEnumerable<Vector<double>> PageRankGenerator(ArrayList at, Vector<double> numLinks, ArrayList leafNodes, double alpha, double convergence, int checkSteps)
	{
		int N = at.Count;
		int M = leafNodes.Count;

		Vector<double> iNew = Ones(N) / N;
		Vector<double> iOld = Ones(N) / N;

		bool done = false;

		while (!done)
		{
			// normalize every now and then for numerical stability
			iNew /= Sum(iNew);

			for (int i = 0; i < checkSteps; i++)
			{
				// swap arrays
				Vector<double> temp = iOld;
				iOld = iNew;
				iNew = temp;

				// an element in the 1 x I vector. 
				// all elements are identical.
				double oneIv = (1 - alpha) * Sum(iOld) / N;

				// an element of the A x I vector.
				// all elements are identical.
				double oneAv = 0.0;
				if (M > 0)
					oneAv = alpha * Sum(Take(iOld, leafNodes)) / N;

				// the elements of the H x I multiplication
				for(int j = 0;j < N;j++)
				{
					List<int> page = (List<int>)at[j];
					double h = 0;

					if (page.Count > 0)
						h = alpha * Take(iOld, page).DotProduct(1.0 / Take(numLinks, page));

					iNew[j] = h + oneAv + oneIv;
				}
			}
			Vector<double> diff = iNew - iOld;
			done = diff.SumMagnitudes() < convergence;

			yield return iNew;
		}
	}

	private Vector<double> Ones(int n)
	{
		Vector<double> result = new DenseVector(n);
		for (int i = 0; i < result.Count; i++)
			result[i] = 1.0;

		return result;
	}

	private double Sum(Vector<double> vector)
	{
		double sum = 0;
		foreach (double val in vector)
			sum += val;

		return sum;
	}

   
	/// <summary>
	/// Simplified (numPy) take method: 1) axis is always 0, 2) first argument is always a vector
	/// </summary>
	/// <param name="vector1">List of values</param>
	/// <param name="vector2">List of indices</param>
	/// <returns>Vector containing elements from vector 1 at the indicies in vector 2</returns>
	private Vector<double> Take(Vector<double> vector1, IList vector2)
	{
		Vector<double> result = new DenseVector(vector2.Count);
		for (int i = 0; i < vector2.Count; i++)
			result[i] = vector1[Convert.ToInt32(vector2[i])];

		return result;
	}

	#endregion
}