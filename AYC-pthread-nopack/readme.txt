                                            ___,,___
                                        ,d8888888888b,_
                                    _,d889'        8888b,
                                _,d8888'          8888888b,
                            _,d8889'           888888888888b,_
                        _,d8889'             888888889'688888, /b
                    _,d8889'               88888889'     `6888d 6,_
                 ,d88886'              _d888889'           ,8d  b888b,  d\
               ,d889'888,             d8889'               8d   9888888Y  )
             ,d889'   `88,          ,d88'                 d8    `,88aa88 9
            d889'      `88,        ,88'                   `8b     )88a88'
           d88'         `88       ,88                   88 `8b,_ d888888
          d89            88,      88                  d888b  `88`_  8888
          88             88b      88                 d888888 8: (6`) 88')
          88             8888b,   88                d888aaa8888, `   'Y'
          88b          ,888888888888                 `d88aa `88888b ,d8
          `88b       ,88886 `88888888                 d88a  d8a88` `8/
           `q8b    ,88'`888  `888'"`88          d8b  d8888,` 88/ 9)_6
             88  ,88"   `88  88p    `88        d88888888888bd8( Z~/
             88b 8p      88 68'      `88      88888888' `688889`
             `88 8        `8 8,       `88    888 `8888,   `qp'
               8 8,        `q 8b       `88  88"    `888b
               q8 8b        "888        `8888'
                "888                     `q88b
                                          "888'

====================================[TEAM]=====================================

                                TEAM ZUBROWKA
                    
                                MEMBERS: 
                                    Jakob Mund
                                    Sebastian Eder
                     
                                AFFILIATION:
                                    TU Munich
                                    
==================================[PROJECT]====================================

                 _     _ _____  ______ _     _        __   __
                 |_____|   |   |  ____ |_____| |        \_/  
                 |     | __|__ |_____| |     | |_____    |   
         _______ _______ _______        _______ ______         _______
         |______ |       |_____| |      |_____| |_____] |      |______
         ______| |_____  |     | |_____ |     | |_____] |_____ |______
          _____  _______  ______ _______               _______       
         |_____] |_____| |_____/ |_____| |      |      |______ |     
         |       |     | |    \_ |     | |_____ |_____ |______ |_____
                  ______ _______ __   _  _____  _______ _______
                 |  ____ |______ | \  | |     | |  |  | |______
                 |_____| |______ |  \_| |_____| |  |  | |______
 _______  _____  _______  _____  _______  ______ _____ _______  _____  __   _
 |       |     | |  |  | |_____] |_____| |_____/   |   |______ |     | | \  |
 |_____  |_____| |  |  | |       |     | |    \_ __|__ ______| |_____| |  \_|


             [Please start counting words below. (should be 298)]
=================================[THE PROBLEM]=================================

Sequences are strings composed of 4 characters. Compare a reference sequence 
with a number of input sequences. Find substrings of at least a given length, 
that are common between the reference sequence and an input sequence.

===================================[THE IDEA]==================================

Design a highly parallelizable algorithm, using a minimal amount of comparisons

================================[THE ALGORITHM]================================

0.  Read input files.

1.  Split input sequences into chunks smaller than the minimum matching length.
    Valid matching sequences must therefore cross at least one chunk's border.
    
2.  Compare only small substrings at the chunks' borders.

3.  Discard chunks that cannot contain matching sequences of the minimum size.
    Having big minimum matching sizes therefore leads to great performance 
    improvements, because large parts of the input sequences are skipped.

4.  From found matches, scan sequences for longer matches. Proceeding from
    the last found match(+1) ensures no unnecessary comparisons are made.

5.  Remove sequences that are enclosed in other sequences.

6.  Order results.

======================[PARALLELIZATION AND OPTIMIZATION]========================

Input: Read input files using a thread pool based on fast spin-locks and 
    paritioned, shared arrays. Employ SSE PCMPISTRI instructions for fast 
    parsing (task-level parallelism, SIMD).

Comparisons: For comparisons, split the diagonals in a number of jobs, where 
    this number is a constant factor wrt the number of threads. Thus, constant 
    synchronization overhead is guaranteed. Furthermore, use SSE 4.2 STTNI 
    instructions for comparisons using inline-assembly (fast xmm0+ register 
    accesses, task-level parallelism).
    
Tasks: Scanning for matches can be a long or short task, depending on the 
    matches' locality. Longer jobs are queued first, then smaller jobs are 
    assigned to later threads.

Removal: The removal of sequences was decoupled to remove data 
    dependencies. The removal is conducted in 2 stages: Scan 
    items and mark for removal (perfectly parallelizable); remove them in 
    linear time (one thread per sequence).

================================================================================

