import random, re, os
import pandas as pd
import collections

with open('../inputs/FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)


def find_paragraphs(txt_directory, terms):
    res = {}
    terms = list(map((lambda x: x.lower()), terms))
    for path, dirs, files in os.walk(txt_directory):
        for file in files:
            full_file = '/'.join([path, file])
            with open(full_file, 'r') as f:
                text = f.read()
            paragraphs = re.split(r'[ \t\r\f\v]*\n[ \t\r\f\v]*\n[ \t\r\f\v]*', text)
            paragraphs = [paragraph.lower() for paragraph in paragraphs if isinstance(paragraph, str)]
            key_paragraphs = []
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs += [paragraph for paragraph in paragraphs
                                   if re.search(re_term, paragraph) and paragraph not in key_paragraphs]
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs = [str(re.sub(re_term, '@@@@' + term, paragraph)) for paragraph in key_paragraphs]
            res[int(file.replace('.txt', ''))] = key_paragraphs
    return res

fcp = ['fcon_1000.projects.nitrc.org', 'Rockland Sample', '1000 Functional Connectomes',
       'International Neuroimaging Data-Sharing Initiative', 'Autism Brain Imaging Data Exchange', 'ADHD-200',
       'Consortium for Reproducibility and Reliability', 'FCP', 'ADHD 200', 'FCON 1000',
       'Functional Connectomes Project', 'www.nitrc.org/projects/fcon_1000', 'NITRC']
paragraph_dict = find_paragraphs('../outputs/txts', fcp)
assign = list(data['i'])
assign_copy = list(assign)


class Member(object):
    def __init__(self, name):
        self.name = name
        self.path = '../outputs/Assignments/' + name + '.txt'
        if os.path.exists(self.path):
            file = open(self.path)
            for i, line in enumerate(file):
                if i == 3:
                    self.articles = [int(l) for l in str(line).strip().split(',')]
                    break
        else:
            self.articles = []

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in self.articles]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.iloc[key]['Title']) +
                        '\n' + str(data.iloc[key]['Authors']) + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraph_dict.items()) if key in self.articles])

    def assign(self, assign, length):
        while length > 0:
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1


class Assignment(object):
    def __init__(self, members):
        self.members = [Member(member) for member in members]

    def assign(self, assignment, length=None, duplicate=True):
        if length is None:
            length = len(assignment)
        length /= len(self.members)
        for member in self.members:
            if not duplicate:
                if len(member.articles) > length:
                    continue
            member.assign(assignment, length)

    def write(self):
        for member in self.members:
            with open(member.path, 'w') as f:
                f.write(str(member))

    def test(self):
        all_articles = []
        for member in self.members:
            all_articles += member.articles
        test_list = [item for item, count in collections.Counter(all_articles).items() if count > 1]
        print(len(test_list))
        print(test_list)


michael = Member('michael')

m_list = """
2
16
20
23
26
29
32
49
60
63
65
69
70
75
80
94
99
103
104
110
115
123
129
130
132
152
155
157
158
161
162
167
172
175
177
184
185
187
190
191
192
194
195
200
205
206
210
217
221
229
232
236
240
244
246
248
250
251
253
261
270
274
275
277
279
283
284
288
290
291
296
299
306
308
311
315
324
330
335
342
353
354
355
367
370
372
383
384
387
395
396
403
416
417
418
423
426
428
430
439
441
442
446
449
452
462
464
465
474
476
481
486
488
491
493
494
499
503
505
507
519
525
526
544
545
550
553
554
560
562
565
577
581
583
589
590
598
600
602
604
605
606
607
609
612
614
621
624
627
629
630
632
634
639
641
655
656
657
666
667
668
669
671
673
674
676
682
683
690
700
701
703
706
710
712
714
715
721
723
725
739
741
744
750
752
753
755
756
760
763
764
767
769
770
772
775
777
779
789
790
798
803
804
805
807
809
812
818
821
823
828
829
835
841
846
849
850
851
854
855
857
858
861
866
868
874
876
878
883
884
891
893
899
900
901
902
904
908
910
916
919
923
929
934
945
946
947
949
963
964
966
970
981
988
990
1001
1002
1005
1006
1007
1012
1015
1019
1022
1029
1032
1043
1044
1045
1048
1054
1055
1059
1064
1069
1070
1073
1082
1086
1095
1096
1099
1104
1105
1111
1114
1116
1132
1143
1148
1150
1157
1164
1167
1169
1175
1177
1192
1195
1201
1204
1205
1206
1210
1212
1213
1214
1223
1226
1230
1240
1241
1253
1258
1267
1268
1269
1271
1274
1280
1281
1290
1292
1298
1303
1310
1312
1313
1316
1320
1325
1329
1331
1335
1336
1337
1342
1345
1359
1365
1367
1376
1379
1381
1382
1391
1394
1395
1400
1403
1407
1411
1419
1423
1430
1435
1437
1441
1442
1446
1447
1451
1455
1456
1459
1462
1465
1466
1470
1474
1485
1489
1491
1495
1496
1497
1498
1512
1521
1522
1523
1525
1531
1539
1552
1554
"""

m_list = [int(i) for i in m_list.strip().split('\n')]

print([i for i in michael.articles if i not in m_list])