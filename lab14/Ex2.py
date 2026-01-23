import math
import string
import matplotlib.pyplot as plt

text_eminescu_train = """
Somnoroase păsărele
Pe la cuiburi se adună,
Se ascund în rămurele -
Noapte bună!

Doar izvoarele suspină,
Pe când codrul negru tace;
Dorm și florile-n grădină -
Dormi în pace!

Trece lebăda pe ape
Între trestii să se culce -
Fie-ți îngerii aproape,
Somnul dulce!

Peste-a nopții feerie
Se ridică mândra lună,
Totu-i vis și armonie -
Noapte bună!
"""

text_stanescu_train = """
Leoaica tânără, iubirea
mi-ai sărit în față.
Mă pândise-n încordare
mai demult.
Colții albi mi i-a înfipt în față,
m-a mușcat leoaica, azi, de față.

Și deodata-n jurul meu, natura
se făcu un cerc, de-a-dura,
când mai larg, când mai aproape,
ca o strângere de ape.
Și privirea-n sus țâșni,
curcubeu tăiat în două,
și auzul o-ntâlni
tocmai lângă ciocârlii.

Mi-am dus mâna la sprânceană,
la tâmplă și la bărbie,
dar mâna nu le mai știe.
Și alunecă în neștire
pe-un deșert în strălucire,
peste care trece-alene
o leoaică arămie
cu mișcările viclene,
încă-o vreme,
și-ncă-o vreme...
"""

text_mixed_test = """
Dintr-un bolovan coboară
pasul tău de domnișoară.
Dintr-o frunză verde, pală
pasul tău de domnișoară.
Dintr-o înserare-n seară
pasul tău de domnișoară.

A fost odată ca-n povești,
a fost ca niciodată,
din rude mari împărătești,
o prea frumoasă fată.
Și era una la părinți
și mândră-n toate cele,
cum e Fecioara între sfinți
și luna între stele.

Pasul tău de domnișoară,
dintr-o pasăre amară.
"""

def preprocess(text):
    text = text.lower()
    text = text.replace('\n', ' ')
    table = str.maketrans('', '', string.punctuation)
    text = text.translate(table)
    return text.split()

words_eminescu = preprocess(text_eminescu_train)
words_stanescu = preprocess(text_stanescu_train)
words_test = preprocess(text_mixed_test)

vocabulary = sorted(list(set(words_eminescu).union(set(words_stanescu)).union(set(words_test))))
word_to_id = {word: i for i, word in enumerate(vocabulary)}
id_to_word = {i: word for word, i in word_to_id.items()}

ids_eminescu = [word_to_id[w] for w in words_eminescu]
ids_stanescu = [word_to_id[w] for w in words_stanescu]
ids_test = [word_to_id[w] for w in words_test]

def train_markov_model(ids, vocab_size):
    counts = {i: {j: 1.0 for j in range(vocab_size)} for i in range(vocab_size)}
    
    for k in range(len(ids) - 1):
        curr_id = ids[k]
        next_id = ids[k+1]
        counts[curr_id][next_id] += 1
            
    probs = {i: {} for i in range(vocab_size)}
    for i in range(vocab_size):
        total = sum(counts[i].values())
        for j in range(vocab_size):
            probs[i][j] = counts[i][j] / total
            
    return probs

vocab_size = len(vocabulary)
prob_eminescu = train_markov_model(ids_eminescu, vocab_size)
prob_stanescu = train_markov_model(ids_stanescu, vocab_size)

def get_log_likelihood_score(id1, id2):
    p_e = prob_eminescu[id1][id2]
    p_s = prob_stanescu[id1][id2]
    return math.log2(p_e / p_s)

window_size = 5
scores = []
x_axis = []

for i in range(len(ids_test) - window_size):
    window = ids_test[i : i + window_size]
    window_score = 0
    
    for j in range(len(window) - 1):
        id_curr = window[j]
        id_next = window[j+1]
        window_score += get_log_likelihood_score(id_curr, id_next)
    
    scores.append(window_score)
    x_axis.append(i)

plt.figure(figsize=(12, 6))
plt.plot(x_axis, scores, label='Log-Likelihood Score', color='purple')
plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.fill_between(x_axis, scores, 0, where=[s > 0 for s in scores], facecolor='blue', alpha=0.3, label='Likely Eminescu (+)')
plt.fill_between(x_axis, scores, 0, where=[s < 0 for s in scores], facecolor='red', alpha=0.3, label='Likely Stanescu (-)')
plt.title("Stylometric Analysis: Eminescu vs. Stanescu (ID Based)")
plt.xlabel("Word Position in Test Text")
plt.ylabel("Log-Likelihood Score")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

print(f"{'Position':<10} {'Window IDs':<25} {'Text Segment':<40} {'Prediction':<10}")
print("-" * 90)

for i in range(0, len(scores), 10): 
    val = scores[i]
    if val > 0.5: detection = "Eminescu"
    elif val < -0.5: detection = "Stanescu"
    else: detection = "Neutral"
    
    window_ids = ids_test[i:i+5]
    snippet_ids = str(window_ids)
    snippet_text = " ".join([id_to_word[idx] for idx in window_ids])
    
    print(f"{i:<10} {snippet_ids:<25} {snippet_text:<40} {detection:<10}")