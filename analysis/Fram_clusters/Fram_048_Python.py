import numpy as np
import pandas as pd

from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.multiclass import OneVsRestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import f1_score, accuracy_score, classification_report, confusion_matrix

# === 1) чтение и базовая подготовка ===
df = pd.read_excel("../../results/Fram/HG_20092023/ASVNoMetazoa/EnvSelected_ClustID.xlsx")
feats = ["doy","Depth","temp","cond","salin","mod_aq_par","mod_aq_par_depth",
         "nitrate_imputed","nitrite_imputed","ammonium_imputed",
         "silicate_imputed","phosphate_imputed"]
target = "ClustID"
group  = "Year"

df = df.dropna(subset=[target]).copy()
df[target] = df[target].astype(int)

# инженерия парочки устойчивых химических индексов + log1p
def add_features(Xdf: pd.DataFrame):
    Z = Xdf.copy()
    eps = 1e-9
    Z["np_ratio"] = Z["nitrate_imputed"] / (Z["phosphate_imputed"] + eps)
    Z["si_n_ratio"] = Z["silicate_imputed"] / (Z["nitrate_imputed"] + eps)
    Z["nh4_nox"] = Z["ammonium_imputed"] / (Z["nitrate_imputed"] + Z["nitrite_imputed"] + eps)
    for c in ["nitrate_imputed","nitrite_imputed","ammonium_imputed","silicate_imputed","phosphate_imputed",
              "Depth","mod_aq_par","mod_aq_par_depth"]:
        Z[c+"_log1p"] = np.log1p(np.clip(Z[c], a_min=0, a_max=None))
    return Z

X_base = df[feats]
X = add_features(X_base)
y = df[target].values
groups = df[group].values

# === 2) модель: OvR-логистическая + калибровка ===
base_clf = OneVsRestClassifier(
    LogisticRegression(
        class_weight="balanced",        # базовая балансировка
        max_iter=3000,
        solver="lbfgs",                 # стабильный для мультинома/ovr
        n_jobs=None
    )
)

# Важный момент: калибруем вероятности (Platt / sigmoid) внутри каждого train-fold
# через внутренний CV -> используем как шаг пайплайна
calibrated = CalibratedClassifierCV(base_clf, cv=3, method="sigmoid")

pipe = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
    ("clf", calibrated),
])

# === 3) leave-one-year-out CV ===
logo = LeaveOneGroupOut()
labels = np.sort(df[target].unique())
K = len(labels)

y_true_all, y_pred_all = [], []
probs_all = []

for tr, te in logo.split(X, y, groups):
    Xtr, Xte = X.iloc[tr], X.iloc[te]
    ytr, yte = y[tr], y[te]

    # Дополнительная «целевая» развертка весов: усилим 6 и 8, ослабим 9 (опционально)
    # Ставь factors по вкусу после первого прогона:
    class_weights = pd.Series(1.0, index=labels)
    class_weights.loc[6] *= 2.0
    class_weights.loc[8] *= 2.5
    class_weights.loc[9] *= 0.7

    # В LogisticRegression мы уже указали class_weight='balanced'.
    # Если хочешь именно эти коэффициенты, можно преобразовать их в sample_weight:
    # -> просто умножим веса классов на balanced-веса
    # Считаем balanced-веса на трейн-фолде:
    vc = pd.Series(ytr).value_counts()
    balanced = {c: (len(ytr)/(K*vc[c])) for c in vc.index}
    sw = np.array([balanced[c] * class_weights[c] for c in ytr])

    # fit с весами (через imputer нужен трюк: fit_step = pipe.named_steps["clf"], но проще так:)
    Xt = pipe.named_steps if False else None  # заглушка: scikit-learn не даёт sample_weight на весь пайплайн
    # поэтому делаем двумя шагами:
    imputer = SimpleImputer(strategy="median")
    scaler  = StandardScaler()

    Xtr_t = scaler.fit_transform(imputer.fit_transform(Xtr))
    Xte_t = scaler.transform(imputer.transform(Xte))

    calibrated = CalibratedClassifierCV(base_clf, cv=3, method="sigmoid")
    calibrated.fit(Xtr_t, ytr, sample_weight=sw)

    y_proba = calibrated.predict_proba(Xte_t)  # shape [n, K]
    probs_all.append(y_proba)
    y_true_all.append(yte)

# === 4) настройка порогов по классам (в т.ч. поднять для класса 9) ===
y_true = np.concatenate(y_true_all)
probs = np.vstack(probs_all)  # согласуем порядок классов
# сопоставление индексов классов:
class_to_idx = {c:i for i, c in enumerate(calibrated.classes_)}

# Подберём пороги по сетке для макро-F1; для 9 стартуем повыше
init_thresholds = np.array([0.5]*K, dtype=float)
if 9 in class_to_idx:
    init_thresholds[class_to_idx[9]] = 0.6   # старт выше для 9

def predict_with_thresholds(P, thr):
    # правило: если все ниже порога — берём argmax; иначе — первый класс, чей proba>=threshold
    # (или можно брать argmax среди прошедших порог)
    passed = (P >= thr[None, :])
    # если ни один не прошёл — argmax
    argmax = P.argmax(axis=1)
    any_passed = passed.any(axis=1)
    # возьмём класс с макс вероятностью среди прошедших
    idx_passed_max = (P * passed).argmax(axis=1)
    yhat = np.where(any_passed, calibrated.classes_[idx_passed_max], calibrated.classes_[argmax])
    return yhat

best_thr = init_thresholds.copy()
best_f1  = -1

grid = np.linspace(0.35, 0.75, 9)  # разумные пороги
for t9 in grid:
    thr = init_thresholds.copy()
    if 9 in class_to_idx:
        thr[class_to_idx[9]] = t9
    # можно аналогично перебирать пороги для 6 и 8
    y_pred = predict_with_thresholds(probs, thr)
    score = f1_score(y_true, y_pred, average="macro")
    if score > best_f1:
        best_f1, best_thr = score, thr.copy()

# финальные метрики
y_pred = predict_with_thresholds(probs, best_thr)
print("Accuracy:", accuracy_score(y_true, y_pred))
print("Macro-F1:", f1_score(y_true, y_pred, average="macro"))
print(pd.DataFrame(classification_report(y_true, y_pred, zero_division=0, output_dict=True)).T)
print("Best thresholds per class:", dict(zip(calibrated.classes_, best_thr)))
print("Confusion matrix:\n", confusion_matrix(y_true, y_pred, labels=calibrated.classes_))
