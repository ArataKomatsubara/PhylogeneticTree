from Bio import SeqIO
from Bio.Phylo import draw_ascii, draw
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO
from Bio import AlignIO
import matplotlib.pyplot as plt
# Fastaファイルを読み込む
alignments= AlignIO.read(r"C:\Users\User\Desktop\komatsubara\kmtD\PythonTest\15F11-Me1.phy", 'phylip')
# 距離行列の計算
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignments)
# 進化系統樹の構築
constructor = DistanceTreeConstructor()
#upgmaを使いたい場合はこっち
tree = constructor.upgma(dm)
#n-j法を使いたい場合はこっち
#tree = constructor.nj(dm)

# Innerを含むノード名を削除するラベル関数
def label_func(node):
    if node.is_terminal():
        return node.name
    else:
        return ""
# 進化系統樹の表示
print(tree)
draw_ascii(tree)
# draw関数で図形を作成
fig, ax = plt.subplots(figsize=(8, 12), dpi=100)
draw(tree, axes=ax, label_func=label_func)
