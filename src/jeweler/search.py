"""Functions for finding optimal codes for coded aperture.

To reduce the search space of the brute force search we need to skip redundant
codes. In order to do this we can take advantage of the invariance of the
Fourier transform to translation. This means that codes which are a translation
of another code are redundant.

The following transformations are not redundant:
    rotation
    flips across any axis

"""

import itertools
import logging
import os

import numpy as np
from sage.all import Necklaces
from tqdm import tqdm

from jeweler.bracelet import bracelet_fc
from jeweler.lyndon import LengthLimitedLyndonWords
from jeweler.io import Archiver

__all__ = [
    'lyndon',
    'bfs',
    'bracelet',
    'necklace',
]

logger = logging.getLogger(__name__)


def _find_codes_prototype(K, L, output_dir, objective_function, density):
    """Search module expects a function with this signature.

    Parameters
    ----------
    K : int
        The minimimum code length inclusive
    L : int
        The maximum code length inclusive
    output_dir : path
        Location to put the output files
    objective_function : function
        A function from jeweler.objective where better scores are larger
    density : float
        The sum of the code divided by the length of the code

    """
    pass


def lyndon(K, L, output_dir, objective_function, density=0.5):
    """Search 1D binary Lyndon words of K <= length <= L for the best codes.

    For the 32-bit space, searching only Lydon words reduces the search space
    to only 3.12% of the full search space.
    """
    logger.info(f"Searching lyndon words of lengths {K}..{L}; "
                f"the objective is '{objective_function.__name__}'.")

    # Stats are L + 1 to avoid repeated subtraction inside search loop
    score_best = np.full(L + 1, -np.inf, dtype=np.float32)
    num_allowed_ones = tuple(int(n * density) for n in range(L + 1))
    num_searched_codes = np.zeros(L + 1, dtype=int)

    # Skip many codes by skipping to the longest one with the desired density.
    w = [0] * (L - num_allowed_ones[-1]) + [1] * num_allowed_ones[-1]

    with Archiver(output_dir=output_dir) as f:
        for code in LengthLimitedLyndonWords(2, L, w):
            if len(code) >= K and np.sum(code) == num_allowed_ones[len(code)]:
                num_searched_codes[len(code)] += 1
                score = objective_function(code)
                if score > score_best[len(code)]:
                    score_best[len(code)] = score
                    f.update(objective_function.__name__,
                             code,
                             score,
                             weight=num_allowed_ones[len(code)])


def bfs(L, density=0.5, batch_size=2**25, filename=None):
    """Find the best binary code of length L using a brute force search.

    Parameters
    ----------
    density: float
        The fraction of indices that are one.
    batch_size: int
        The number of bits to try at once. Limits memory consumption.
    filename: string
        The name of a file to dump the best result at the completion of every
        batch.

    """
    k = int(L * density)  # number of 1s in the code
    code_generator = itertools.combinations(range(L), k)
    num_combinations = (np.math.factorial(L) //
                        (np.math.factorial(k) * np.math.factorial(L - k)))
    # Determine the number of codes to try in one go and make a code generator
    num_codes_in_batch = max(batch_size // L, 1)
    num_batches = np.ceil(num_combinations / num_codes_in_batch).astype(int)

    code_best = None
    score_best = -np.inf
    num_searched_codes = 0

    title = "brute force 1D {:d}-bit code".format(L)
    for i in tqdm(range(0, num_batches), desc=title):
        # Get the next batch of codes
        indices = list(itertools.islice(code_generator, num_codes_in_batch))
        if len(indices) == 0:
            # Quit if there are no more codes to try
            break
        num_searched_codes += len(indices)
        # convert indices to a binary code
        codes = np.zeros([len(indices), L], dtype=np.float32)
        rows, cols = np.meshgrid(range(len(indices)), range(k), indexing='ij')
        codes[rows, indices] = 1
        # compute the DFT of the code
        dft = np.fft.rfft(codes, axis=1)
        # rate all of the DFTs
        score = np.min(np.abs(dft), axis=1) - np.var(dft, axis=1)
        best = np.argmax(score)
        if score[best] > score_best:
            code_best = codes[best]
            score_best = score[best]
            if filename is not None:
                with open(filename, mode='w', encoding='utf-8') as f:
                    print("# Best code from searching {:,d} codes".format(
                        num_searched_codes),
                          file=f)
                    print(code_best.astype(int), file=f, flush=True)
    return code_best


def bracelet(
    K,
    L,
    output_dir,
    objective_function,
    density=0.5,
    batch_size=2048,
):
    """Search bracelets of length L and fixed content for the best binary code.

    Parameters
    ----------
    L : int
        The maximum code length
    output_dir : path
        Location to put the output files
    objective_function : function
        A function from jeweler.objective where better scores are larger
    density : float
        The sum of the code divided by the length of the code
    """
    logger.info(f"Fixed-content bracelets of length {K}..{L}.")
    logger.info(f"the objective is '{objective_function.__name__}'.")
    logger.info(f"code density is {density:g}.")

    with Archiver(output_dir=output_dir) as f:

        for L in range(K, L + 1):

            logger.info(f"Generating bracelets of length {L}.")
            k = int(L * density)  # number of 1s in the code
            codes = bracelet_fc(L, 2, [L - k, k])
            logger.info(f"{len(codes):,d} bracelets discovered.")

            title = f"fixed-content bracelets 1D {L:d}-bit code"
            code_best = None
            score_best = -np.inf
            for _ in tqdm(range(len(codes) // batch_size + 1), desc=title):
                batch = codes[-batch_size:]
                del codes[-batch_size:]
                scores = objective_function(batch)
                best = np.argmax(scores)
                if scores[best] > score_best:
                    score_best = scores[best]
                    code_best = batch[best]

            f.update(
                objective_function.__name__,
                code_best,
                score_best,
                weight=k,
            )


def _necklace_chunk(necklaces, chunksize):
    """Wrap a sage.Necklaces instance in a Python generator.

    Parameters
    ----------
    necklages : sage.Necklaces
        A SageMath necklage class instance.
    chunksize : 2-tuple
        The shape of the desired chunk of necklaces.
    """
    assert chunksize[1] == len(necklaces.first())
    chunk = np.empty(chunksize, dtype='float32')
    i = 0
    for x in necklaces:
        chunk[i] = x
        i += 1
        if i >= chunksize[0]:
            yield chunk - 1
            chunk = np.empty(chunksize, dtype='float32')
            i = 0
    if i < chunksize[0]:
        yield chunk[:i] - 1


def necklace(
    K,
    L,
    output_dir,
    objective_function,
    density=0.5,
    batch_size=2048,
):
    """Search necklaces of length L and fixed content for the best binary code.

    Parameters
    ----------
    L : int
        The maximum code length
    output_dir : path
        Location to put the output files
    objective_function : function
        A function from jeweler.objective where better scores are larger
    density : float
        The sum of the code divided by the length of the code
    """
    logger.info(f"Fixed-content necklaces of length {K}..{L}.")
    logger.info(f"the objective is '{objective_function.__name__}'.")
    logger.info(f"code density is {density:g}.")

    with Archiver(output_dir=output_dir) as f:

        for L in range(K, L + 1):

            logger.info(f"Generating necklaces of length {L}.")
            k = int(L * density)  # number of 1s in the code

            codes = Necklaces([L - k, k])
            ncodes = codes.cardinality()
            chunks = _necklace_chunk(codes, (batch_size, L))
            logger.info(f"{ncodes:,d} necklaces discovered.")

            title = f"fixed-content necklaces 1D {L:d}-bit code"
            code_best = None
            score_best = -np.inf
            for _ in tqdm(range(ncodes // batch_size + 1), desc=title):
                batch = next(chunks)
                scores = objective_function(batch)
                best = np.argmax(scores)
                if scores[best] > score_best:
                    score_best = scores[best]
                    code_best = batch[best].astype('int', copy=False).tolist()

            f.update(
                objective_function.__name__,
                code_best,
                score_best,
                weight=k,
            )
