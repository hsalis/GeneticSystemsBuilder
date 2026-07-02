import sys
import shutil
from . import utils

from rocksdict import Rdict, WriteBatch

class kmerSetDB(object):
    '''
    Plyvel based scalable on disk kmerSetDB
    for k-mer storage.
    '''

    def __init__(self, path, homology, verbose=False):
        '''
        kmerSetDB constructor.
        '''
        try:
            # Assert path is string
            if not isinstance(path, str):
                print('\n[Non-Repetitive Parts Calculator - Background]')
                print(('\n [ERROR]    Path must be a string, not {}'.format(
                    path)))
                print(' [SOLUTION] Try correcting Path\n')
                raise ValueError

            # Assert homology is positive integer
            if not isinstance(homology, int):
                print('\n[Non-Repetitive Parts Calculator - Background]')
                print(('\n [ERROR]    Lmax must be an integer, not {}'.format(type(homology))))
                print(' [SOLUTION] Try correcting Lmax\n')
                raise ValueError
            if homology-1 < 5:
                print('\n[Non-Repetitive Parts Calculator - Background]')
                print(('\n [ERROR]    Lmax must be greater than 4, not \'{}\''.format(homology-1)))
                print(' [SOLUTION] Try correcting Lmax\n')
                raise ValueError

            # Path setup
            self.PATH = path.rstrip('/') + '/'

            # Create/Open RocksDB Rdict object
            self.DB = Rdict(self.PATH)

            # Length setup
            try:
                self.LEN = self.DB['LEN']
            except:
                self.LEN = 0
                self.DB['LEN']= 0

            # K setup
            try:
                self.K = int(self.DB['K'])
            except:
                self.K = int(homology)
                self.DB['K'] = int(homology)

            # Verbosity setup
            self.VERB = bool(verbose)

            # Object ALIVE status
            self.ALIVE = True

        except ValueError as E:
            raise ValueError('Invalid Path or Lmax')
        except Exception as E:
            raise E

    def __repr__(self):
        '''
        User function to return the string representation of
        kmerSetDB.
        '''
        return 'kmerSetDB stored at {} with {} {}-mers'.format(
            self.PATH, self.LEN, self.K)

    def __str__(self):
        '''
        See __repr__
        '''
        return self.__repr__()

    def alivemethod(method):
        '''
        Internal decorator to gate kmerSetDB operation
        once dropped.
        '''
        def wrapper(self, *args, **kwargs):
            if self.ALIVE:
                return method(self, *args, **kwargs)
            else:
                raise RuntimeError('kmerSetDB was closed or dropped')
        return wrapper

    def _verb_action(self, action, index, seq):
        '''
        Internal helper function to print sequence
        updates to kmerSetDB based on verbosity.
        '''
        if self.VERB:
            psl = min(len(seq), 10)
            printed = '  {} Seq {}: {}...'.format(
                action, index, seq[:psl])
            clear_length = len(printed)
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            sys.stdout.flush()

    def _get(self, key):
        '''
        Internal helper to fetch value for given
        key, from kmerSetDB instance.
        '''
        try:
            val = self.DB[key]
            return str(val)
        except:
            return None

    @alivemethod
    def add(self, seq):
        '''
        User function to add seq k-mers to kmerSetDB.
        '''
        try:
            if not seq in ['K', 'LEN']:
                seq = seq.replace('U', 'T')
                for kmer in utils.stream_min_kmers(seq=seq, k=self.K):
                    Exists, key = self.DB.key_may_exist(kmer, fetch=True)
                    if not Exists:
                        self.DB[kmer] = True
                        self.LEN += 1
            self.DB['LEN'] = self.LEN
        except Exception as E:
            raise E

    @alivemethod
    def multiadd(self, seq_list):
        '''
        User function to add seq k-mers for each seq in
        seq_list to kmerSetDB via single transaction.
        '''
        #try:
        pt = False
        WB = WriteBatch()
        index = 0
        for seq in seq_list:
            if not pt and self.VERB:
#                    print('\n[Background Processing]')
                pt = True
            if not seq in ['K', 'LEN']:
                self._verb_action(
                    action='Adding',
                    index=index,
                    seq=seq)
                index += 1
                seq = seq.replace('U', 'T')
                for kmer in utils.stream_min_kmers(seq=seq, k=self.K):
                    Exists, key = self.DB.key_may_exist(kmer, fetch=True)
                    if not Exists:
                        WB.put(kmer, True)
                        self.LEN += 1
        WB.put('LEN', self.LEN)
        self.DB.write(WB)
        #except Exception as E:
        #    raise E

    @alivemethod
    def __contains__(self, seq, rna=False):
        '''
        Python dunder function to check existence of
        any k-mer from seq in kmerSetDB.
        '''
        try:
            if not rna:
                seq = seq.replace('U', 'T')
            for kmer in utils.stream_min_kmers(seq=seq, k=self.K):
                Exists, key = self.DB.key_may_exist(kmer, fetch=True)
                if Exists:
                    return True
            return False
        except Exception as E:
            raise E

    @alivemethod
    def multicheck(self, seq_list):
        '''
        User function to check existence of any k-mer for
        each seq in seq_list in kmerSetDB via single
        transaction.
        '''
        try:
            for seq in seq_list:
                yield self.__contains__(
                    seq=seq,
                    rna=False)
        except Exception as E:
            raise E

    @alivemethod
    def __iter__(self):
        '''
        User fuction to iterate over k-mers stored in
        kmerSetDB.
        '''
        for key, _ in self.DB.items():
            if not key in ['K', 'LEN']:
                yield key

    @alivemethod
    def __len__(self):
        '''
        User function to return the number of keys stored
        in kmerSetDB.
        '''
        return self.LEN

    @alivemethod
    def remove(self, seq):
        '''
        User function to remove all k-mers in seq
        from kmerSetDB.
        '''
        try:
            if not seq in ['K', 'LEN']:
                for kmer in utils.stream_min_kmers(seq=seq, k=self.K):
                    Exists, key = self.DB.key_may_exist(kmer, fetch=True)
                    if Exists:
                        self.DB.delete(kmer)
                        self.LEN -= 1
                self.DB['LEN'] = self.LEN
        except Exception as E:
            raise E

    @alivemethod
    def multiremove(self, seq_list, clear=False):
        '''
        User function to remove all k-mers from each
        seq in seq_list from kmerSetDB via single
        transaction.
        '''
        try:
            pt = False
            WB = WriteBatch()
            index = 0
            for seq in seq_list:
                if not pt and self.VERB:
#                    print('\n[Background Processing]')
                    pt = True
                if not seq in ['K', 'LEN']:
                    self._verb_action(
                        action='Removing',
                        index=index,
                        seq=seq)
                    index += 1
                    for kmer in utils.stream_min_kmers(seq=seq, k=self.K):
                        Exists, key = self.DB.key_may_exist(kmer, fetch=True)
                        if Exists:
                            WB.delete(kmer)
                            self.LEN -= 1
            self.DB.write(WB)
            self.DB['LEN'] = self.LEN
            if self.VERB: print()
        except Exception as E:
            raise E

    def clear(self, Lmax=None):
        '''
        User function to clear all k-mers stored in
        kmerSetDB.
        '''
        self.drop()
        self.DB = Rdict(self.PATH)
        self.LEN = 0
        self.DB['LEN'] = self.LEN
        if Lmax is None:
            self.DB['K'] = self.K
        else:
            if not isinstance(Lmax, int):
                print(('\n [ERROR]    Lmax must be an integer, not {}'.format(
                    type(Lmax))))
                print(' [SOLUTION] Try correcting Lmax\n')
                raise ValueError
            self.K = int(Lmax) + 1
            if Lmax < 5:
                print(('\n [ERROR]    Lmax must be greater than 4, not {}'.format(
                    Lmax)))
                print(' [SOLUTION] Try correcting Lmax\n')
                raise ValueError
            self.DB['K'] = self.K
        self.ALIVE = True

    def close(self):
        '''
        User function to close kmerSetDB.
        '''
        if self.ALIVE:
            del self.DB
            self.DB = None
            self.ALIVE = False
            return True
        return False

    def drop(self):
        '''
        User function to drop kmerSetDB.
        '''
        if self.ALIVE:
            del self.DB
            self.DB = None
            shutil.rmtree(self.PATH)
            self.ALIVE = False
            return True
        return False

def test():

    # required test imports
    import random
    import time

    Lmax = 99

    # returns a random DNA string of length
    def get_DNA(length):
        return ''.join(random.choice('ATGC') for _ in range(length))

    # create corpus
    total_elements = 10000
    dna_strings = sorted(get_DNA(100) for _ in range(total_elements))
    print(('Testing Operations on {} items'.format(total_elements*2)))

    # clearance
    try:
        shutil.rmtree('./testBG/')
    except:
        pass

    # object test
    myBG = kmerSetDB(homology=Lmax, path='./testBG')
    assert isinstance(myBG, kmerSetDB)
    assert repr(myBG) == 'kmerSetDB stored at ./testBG/ with 0 {}-mers'.format(Lmax)

    # add test
    t0 = time.time()
    for i in range(total_elements):
        myBG.add(dna_strings[i])
    print(('add          time = {} sec'.format(time.time()-t0)))
    myBG.clear()

    # # multiadd test
    t0 = time.time()
    myBG.multiadd(dna_strings)
    print(('multiadd     time = {} sec'.format(time.time()-t0)))

    # __contains__ test
    seq = get_DNA(length=100)
    assert seq not in myBG
    t0 = time.time()
    for dna in dna_strings:
        assert dna in myBG
    print(('__contains__ time = {} sec'.format(time.time()-t0)))

    # multicheck test
    t0 = time.time()
    assert all(myBG.multicheck(seq_list=dna_strings))
    print(('multicheck   time = {} sec'.format(time.time()-t0)))

    # __iter__ test
    t0 = time.time()
    sorted(myBG)
    print(('__iter__     time = {} sec'.format(time.time()-t0)))

    # __len__ test
    t0 = time.time()
    assert len(myBG) == total_elements * 2
    print(('__len__      time = {} sec'.format(time.time()-t0)))

    # remove test
    t0 = time.time()
    assert len(myBG) == total_elements*2
    for dna_string in dna_strings:
        myBG.remove(dna_string)
    print(('remove       time = {} sec'.format(time.time()-t0)))
    assert len(myBG) == 0

    # multiremove test
    myBG.multiadd(dna_strings)
    assert len(myBG) == total_elements*2
    t0 = time.time()
    myBG.multiremove(dna_strings)
    print(('multiremove  time = {} sec'.format(time.time()-t0)))
    assert len(myBG) == 0

    # clear test
    myBG.multiadd(dna_strings)
    t0 = time.time()
    myBG.clear()
    print(('clear        time = {} sec'.format(time.time()-t0)))

    # close test
    t0 = time.time()
    myBG.close()
    print(('close        time = {} sec'.format(time.time()-t0)))

    # drop test
    myBG = kmerSetDB(homology=Lmax, path='./testBG')
    t0 = time.time()
    myBG.drop()
    print(('drop         time = {} sec'.format(time.time()-t0)))

if __name__ == '__main__':
    test()