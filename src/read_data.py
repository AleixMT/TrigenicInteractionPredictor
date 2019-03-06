import sys

info = open(sys.argv[1]).readlines()
percent = 0.2
line = 0
predictions = []
triplets = {}
positives = []
negatives = []
count = 0
maxL = 0
for line in info:
    data = line.strip().split('\t')

    if data[0]=='Max Likelihood:': maxL=float(data[1])
    if data==[''] and count == 1: break
    if count == 1:
        score,triplet,interaction = data[0],data[1],data[2]
        predictions.append(triplet)
        triplets[triplet] = score
        if int(interaction)==1:
            positives.append(triplet)
        else:
            negatives.append(triplet)
    if data[0]=='Predicted Interaction': count = 1


#make two vectors of ones and zeroes for precision and recall
nedges1 = int(len(predictions)*percent)
ones = predictions[:nedges1]
zeroes = predictions[nedges1 : ]


##compute metrics

##auc
auc = 0
for pos in positives:
    for neg in negatives:
        if float(triplets[pos]) > float(triplets[neg]):
            auc += 1
        elif float(triplets[pos]) == float(triplets[neg]):
            auc += .5
if auc > 0: auc/=float(len(positives)*len(negatives))

truepositives = 0
precision = 0
recall = 0
for one in ones:
    if one in positives:
        truepositives += 1
        
if truepositives > 0:
    precision=truepositives/float(len(ones))
    recall=truepositives/float(len(positives))


print sys.argv[1],maxL,auc,precision,recall,len(positives)/float(len(positives)+len(negatives)) 
