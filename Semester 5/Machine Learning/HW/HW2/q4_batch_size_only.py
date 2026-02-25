"""
Q4: Batch Size Experiment Only

This script runs only the batch size comparison experiment from the notebook.
Uses lr=5e-4, filter_size=5, num_conv_blocks=4 (matching the notebook configuration).
"""

import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision.transforms as transforms
import torchvision.transforms.functional as TF
from torch.utils.data import DataLoader
from datasets import load_dataset
import pytorch_lightning as pl
import torchmetrics
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Set random seeds
torch.manual_seed(67)
np.random.seed(67)
random.seed(67)

# Device configuration
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device}")
torch.set_float32_matmul_precision('high')

# Dataset + simple transform (with light augmentation for train)
IMG_SIZE = 224
MEAN = [0.485, 0.456, 0.406]
STD = [0.229, 0.224, 0.225]

_color_jitter = transforms.ColorJitter(brightness=0.15, contrast=0.15, saturation=0.15)

def _scale_bbox(box, orig_w, orig_h, new_w, new_h):
    x1, y1, x2, y2 = box
    sx = new_w / orig_w
    sy = new_h / orig_h
    return [x1 * sx, y1 * sy, x2 * sx, y2 * sy]

def _transform_one(img, bbox, train=False):
    orig_w, orig_h = img.size
    img = TF.resize(img, (IMG_SIZE, IMG_SIZE))
    bbox = _scale_bbox(bbox, orig_w, orig_h, IMG_SIZE, IMG_SIZE)

    if train and random.random() < 0.5:
        img = TF.hflip(img)
        x1, y1, x2, y2 = bbox
        bbox = [IMG_SIZE - x2, y1, IMG_SIZE - x1, y2]

    if train and random.random() < 0.3:
        img = _color_jitter(img)

    img = TF.to_tensor(img)
    img = TF.normalize(img, mean=MEAN, std=STD)
    bbox = torch.tensor([b / IMG_SIZE for b in bbox], dtype=torch.float32)
    return img, bbox

def transform_example(example, train=False):
    imgs = example["img"]
    bboxes = example["bbox"]
    if isinstance(imgs, list):
        new_imgs = []
        new_boxes = []
        for i, img in enumerate(imgs):
            img_t, box_t = _transform_one(img, bboxes[i], train=train)
            new_imgs.append(img_t)
            new_boxes.append(box_t)
        example["img"] = new_imgs
        example["bbox"] = new_boxes
    else:
        img_t, box_t = _transform_one(imgs, bboxes, train=train)
        example["img"] = img_t
        example["bbox"] = box_t
    return example

def transform_train(example):
    return transform_example(example, train=True)

def transform_eval(example):
    return transform_example(example, train=False)

_ds = load_dataset("cvdl/oxford-pets")

# Keep only needed columns to avoid collate issues
keep_cols = {"img", "bbox", "category"}
for split_name in ["train", "valid", "test"]:
    cols = list(_ds[split_name].column_names)
    drop_cols = [c for c in cols if c not in keep_cols]
    if drop_cols:
        _ds[split_name] = _ds[split_name].remove_columns(drop_cols)

train_ds = _ds["train"].with_transform(transform_train)
valid_ds = _ds["valid"].with_transform(transform_eval)
test_ds = _ds["test"].with_transform(transform_eval)

print("Dataset loaded")

# IoU and mAP metric functions
def compute_iou(pred_bbox, gt_bbox):
    """Compute Intersection over Union between predicted and ground truth bounding boxes."""
    pred_x1, pred_y1, pred_x2, pred_y2 = pred_bbox
    gt_x1, gt_y1, gt_x2, gt_y2 = gt_bbox
    
    pred_x1, pred_x2 = min(pred_x1, pred_x2), max(pred_x1, pred_x2)
    pred_y1, pred_y2 = min(pred_y1, pred_y2), max(pred_y1, pred_y2)
    gt_x1, gt_x2 = min(gt_x1, gt_x2), max(gt_x1, gt_x2)
    gt_y1, gt_y2 = min(gt_y1, gt_y2), max(gt_y1, gt_y2)
    
    inter_x1 = max(pred_x1, gt_x1)
    inter_y1 = max(pred_y1, gt_y1)
    inter_x2 = min(pred_x2, gt_x2)
    inter_y2 = min(pred_y2, gt_y2)
    
    inter_area = max(0, inter_x2 - inter_x1) * max(0, inter_y2 - inter_y1)
    
    pred_area = (pred_x2 - pred_x1) * (pred_y2 - pred_y1)
    gt_area = (gt_x2 - gt_x1) * (gt_y2 - gt_y1)
    union_area = pred_area + gt_area - inter_area
    
    iou = inter_area / union_area if union_area > 0 else 0
    return iou

def compute_batch_iou(pred_bboxes, gt_bboxes):
    """Compute IoU for a batch of bounding boxes."""
    ious = []
    for i in range(pred_bboxes.shape[0]):
        pred = pred_bboxes[i].cpu().numpy() if torch.is_tensor(pred_bboxes[i]) else pred_bboxes[i]
        gt = gt_bboxes[i].cpu().numpy() if torch.is_tensor(gt_bboxes[i]) else gt_bboxes[i]
        iou = compute_iou(pred, gt)
        ious.append(iou)
    return torch.tensor(ious)

def compute_ap(precisions, recalls):
    """Compute Average Precision using 11-point interpolation."""
    ap = 0.0
    for t in np.arange(0.0, 1.1, 0.1):
        if np.sum(recalls >= t) == 0:
            p = 0
        else:
            p = np.max(precisions[recalls >= t])
        ap += p / 11.0
    return ap

def compute_map(pred_labels, pred_bboxes, gt_labels, gt_bboxes, num_classes, iou_threshold=0.5):
    """Compute mean Average Precision for object detection."""
    ious = compute_batch_iou(pred_bboxes, gt_bboxes)
    
    aps = []
    for cls in range(num_classes):
        gt_mask = gt_labels == cls
        pred_mask = pred_labels == cls
        
        if not gt_mask.any():
            continue
        
        tp = (pred_mask & gt_mask & (ious > iou_threshold)).float()
        fp = (pred_mask & (~gt_mask | (ious <= iou_threshold))).float()
        
        tp_cumsum = torch.cumsum(tp, dim=0)
        fp_cumsum = torch.cumsum(fp, dim=0)
        
        recalls = tp_cumsum / gt_mask.sum().float()
        precisions = tp_cumsum / (tp_cumsum + fp_cumsum + 1e-10)
        
        ap = compute_ap(precisions.cpu().numpy(), recalls.cpu().numpy())
        aps.append(ap)
    
    return np.mean(aps) if aps else 0.0

print("Metrics functions ready")

# FlexibleCNN Model
class FlexibleCNN(pl.LightningModule):
    def __init__(self, num_classes=37, lr=1e-3, bbox_loss_weight=3.0, bbox_lr_mult=2.0, 
                 num_conv_blocks=4, filter_size=3, base_channels=32):
        super().__init__()
        self.save_hyperparameters()
        self.lr = lr
        self.bbox_loss_weight = bbox_loss_weight
        self.bbox_lr_mult = bbox_lr_mult
        self.num_classes = num_classes
        
        layers = []
        in_channels = 3
        out_channels = base_channels
        padding = filter_size // 2
        
        for i in range(num_conv_blocks):
            layers.extend([
                nn.Conv2d(in_channels, out_channels, kernel_size=filter_size, padding=padding),
                nn.BatchNorm2d(out_channels),
                nn.ReLU(),
                nn.Conv2d(out_channels, out_channels, kernel_size=filter_size, padding=padding),
                nn.BatchNorm2d(out_channels),
                nn.ReLU(),
                nn.MaxPool2d(2),
            ])
            in_channels = out_channels
            out_channels = min(out_channels * 2, 512)
        
        self.features = nn.Sequential(*layers)
        final_channels = in_channels
        
        self.pool = nn.AdaptiveAvgPool2d((4, 4))
        flat_size = final_channels * 4 * 4
        
        self.dropout = nn.Dropout(0.4)
        self.class_head = nn.Sequential(
            nn.Linear(flat_size, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, num_classes)
        )
        self.bbox_head = nn.Sequential(
            nn.Linear(flat_size, 512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 4)
        )
        
        self.train_acc = torchmetrics.Accuracy(task="multiclass", num_classes=num_classes)
        self.val_acc = torchmetrics.Accuracy(task="multiclass", num_classes=num_classes)
        self.test_acc = torchmetrics.Accuracy(task="multiclass", num_classes=num_classes)
        
        self.val_ious = []
        self.val_pred_labels = []
        self.val_pred_bboxes = []
        self.val_gt_labels = []
        self.val_gt_bboxes = []
        
        self.test_ious = []
        self.test_pred_labels = []
        self.test_pred_bboxes = []
        self.test_gt_labels = []
        self.test_gt_bboxes = []
        
        # Track epoch-level losses for visualization
        self.train_losses = []
        self.val_losses = []

    def forward(self, x):
        feats = self.features(x)
        pooled = self.pool(feats).view(x.size(0), -1)
        class_feats = self.dropout(pooled)
        class_logits = self.class_head(class_feats)
        bbox = self.bbox_head(pooled)
        return class_logits, bbox
    
    def _unpack(self, batch):
        if isinstance(batch, dict):
            images = batch["img"]
            labels = batch["category"]
            bboxes = batch["bbox"]
        else:
            images, labels, bboxes = batch
        if not torch.is_tensor(bboxes):
            bboxes = torch.tensor(bboxes, device=images.device, dtype=torch.float32)
        bboxes = bboxes.float()
        labels = labels.long()
        return images, labels, bboxes

    def _shared_step(self, batch):
        images, labels, bboxes = self._unpack(batch)
        logits, pred_bboxes = self(images)
        loss_cls = F.cross_entropy(logits, labels)
        pred_bboxes_norm = torch.sigmoid(pred_bboxes)
        loss_bbox = F.smooth_l1_loss(pred_bboxes_norm, bboxes)
        loss = loss_cls + self.bbox_loss_weight * loss_bbox
        preds = torch.argmax(logits, dim=1)
        return loss, loss_cls, loss_bbox, preds, labels, pred_bboxes_norm, bboxes

    def training_step(self, batch, batch_idx):
        loss, loss_cls, loss_bbox, preds, labels, pred_bboxes_norm, bboxes = self._shared_step(batch)
        self.train_acc(preds, labels)
        self.log("train_loss", loss, prog_bar=True)
        self.log("train_acc", self.train_acc, prog_bar=True)
        self.log("train_cls_loss", loss_cls)
        self.log("train_bbox_loss", loss_bbox)
        return loss
    
    def on_train_epoch_end(self):
        # Store epoch-level training loss for plotting
        if 'train_loss' in self.trainer.logged_metrics:
            self.train_losses.append(self.trainer.logged_metrics['train_loss'].item())

    def validation_step(self, batch, batch_idx):
        loss, loss_cls, loss_bbox, preds, labels, pred_bboxes_norm, bboxes = self._shared_step(batch)
        self.val_acc(preds, labels)
        ious = compute_batch_iou(pred_bboxes_norm, bboxes)
        mean_iou = ious.mean()
        self.val_ious.extend(ious.tolist())
        self.val_pred_labels.append(preds.cpu())
        self.val_pred_bboxes.append(pred_bboxes_norm.cpu())
        self.val_gt_labels.append(labels.cpu())
        self.val_gt_bboxes.append(bboxes.cpu())
        self.log("val_loss", loss, prog_bar=True)
        self.log("val_acc", self.val_acc, prog_bar=True)
        self.log("val_cls_loss", loss_cls)
        self.log("val_bbox_loss", loss_bbox)
        self.log("val_iou", mean_iou, prog_bar=True)

    def on_validation_epoch_end(self):
        if len(self.val_pred_labels) > 0:
            pred_labels = torch.cat(self.val_pred_labels)
            pred_bboxes = torch.cat(self.val_pred_bboxes)
            gt_labels = torch.cat(self.val_gt_labels)
            gt_bboxes = torch.cat(self.val_gt_bboxes)
            map_score = compute_map(pred_labels, pred_bboxes, gt_labels, gt_bboxes, 
                                   self.num_classes, iou_threshold=0.5)
            mean_iou = np.mean(self.val_ious)
            self.log("val_map", map_score, prog_bar=True)
            self.log("val_mean_iou", mean_iou)
            
            # Store epoch-level validation loss for plotting
            if 'val_loss' in self.trainer.logged_metrics:
                self.val_losses.append(self.trainer.logged_metrics['val_loss'].item())
            
            self.val_ious = []
            self.val_pred_labels = []
            self.val_pred_bboxes = []
            self.val_gt_labels = []
            self.val_gt_bboxes = []

    def test_step(self, batch, batch_idx):
        loss, loss_cls, loss_bbox, preds, labels, pred_bboxes_norm, bboxes = self._shared_step(batch)
        self.test_acc(preds, labels)
        ious = compute_batch_iou(pred_bboxes_norm, bboxes)
        mean_iou = ious.mean()
        self.test_ious.extend(ious.tolist())
        self.test_pred_labels.append(preds.cpu())
        self.test_pred_bboxes.append(pred_bboxes_norm.cpu())
        self.test_gt_labels.append(labels.cpu())
        self.test_gt_bboxes.append(bboxes.cpu())
        self.log("test_loss", loss, prog_bar=True)
        self.log("test_acc", self.test_acc, prog_bar=True)
        self.log("test_cls_loss", loss_cls)
        self.log("test_bbox_loss", loss_bbox)
        self.log("test_iou", mean_iou, prog_bar=True)

    def on_test_epoch_end(self):
        if len(self.test_pred_labels) > 0:
            pred_labels = torch.cat(self.test_pred_labels)
            pred_bboxes = torch.cat(self.test_pred_bboxes)
            gt_labels = torch.cat(self.test_gt_labels)
            gt_bboxes = torch.cat(self.test_gt_bboxes)
            map_score = compute_map(pred_labels, pred_bboxes, gt_labels, gt_bboxes, 
                                   self.num_classes, iou_threshold=0.5)
            mean_iou = np.mean(self.test_ious)
            self.log("test_map", map_score, prog_bar=True)
            self.log("test_mean_iou", mean_iou)
            self.test_ious = []
            self.test_pred_labels = []
            self.test_pred_bboxes = []
            self.test_gt_labels = []
            self.test_gt_bboxes = []

    def configure_optimizers(self):
        return torch.optim.Adam(
            [
                {"params": self.features.parameters(), "lr": self.lr},
                {"params": self.class_head.parameters(), "lr": self.lr},
                {"params": self.bbox_head.parameters(), "lr": self.lr * self.bbox_lr_mult},
            ]
        )

print("FlexibleCNN model defined")


if __name__ == "__main__":
    # Configuration matching the notebook
    batch_sizes = [16, 32, 64, 128]
    NUM_WORKERS = 4  # Enable multiple workers for speed
    batch_results = []
    training_histories = {}  # Store training history for loss curves

    print("\n" + "="*70)
    print("BATCH SIZE EXPERIMENT")
    print("="*70)
    print(f"Configuration: lr=5e-4, filter_size=5, num_conv_blocks=4")
    print(f"Batch sizes: {batch_sizes}")
    print(f"Workers: {NUM_WORKERS}")
    print("="*70 + "\n")
    
    for batch_size in batch_sizes:
        print(f"\n{'='*50}")
        print(f"Training with batch size: {batch_size}")
        print('='*50)
        
        # Create new data loaders with different batch size
        train_loader_exp = DataLoader(train_ds, batch_size=batch_size, shuffle=True, 
                                      num_workers=NUM_WORKERS, pin_memory=torch.cuda.is_available())
        valid_loader_exp = DataLoader(valid_ds, batch_size=batch_size, shuffle=False, 
                                      num_workers=NUM_WORKERS, pin_memory=torch.cuda.is_available())
        test_loader_exp = DataLoader(test_ds, batch_size=batch_size, shuffle=False, 
                                     num_workers=NUM_WORKERS, pin_memory=torch.cuda.is_available())
        
        # Model configuration matching the notebook
        model = FlexibleCNN(
            num_classes=37, 
            lr=5e-4,  # Matching notebook
            bbox_loss_weight=5.0, 
            bbox_lr_mult=2,
            num_conv_blocks=4,
            filter_size=5  # Matching notebook
        )
        
        trainer = pl.Trainer(
            max_epochs=50,
            accelerator=device,
            log_every_n_steps=10,
            enable_progress_bar=True,
            enable_model_summary=False,
            default_root_dir=f"lightning_logs/batch_size_{batch_size}"
        )
        
        trainer.fit(model, train_loader_exp, valid_loader_exp)
        test_results = trainer.test(model, test_loader_exp, verbose=False)
        
        # Extract training history
        train_losses = []
        val_losses = []
        for metric_dict in trainer.logged_metrics:
            if 'train_loss' in metric_dict:
                train_losses.append(metric_dict['train_loss'])
            if 'val_loss' in metric_dict:
                val_losses.append(metric_dict['val_loss'])
        
        # Store history with unique key
        training_histories[batch_size] = {
            'train_loss': train_losses if train_losses else [],
            'val_loss': val_losses if val_losses else [],
            'model': model
        }
        
        result = {
            'batch_size': batch_size,
            'test_acc': test_results[0]['test_acc'],
            'test_iou': test_results[0].get('test_mean_iou', 0),
            'test_map': test_results[0].get('test_map', 0),
            'test_loss': test_results[0]['test_loss']
        }
        batch_results.append(result)
        print(f"\nResults: Acc={result['test_acc']:.4f}, IoU={result['test_iou']:.4f}, mAP={result['test_map']:.4f}")

    # Visualize results
    print("\n" + "="*70)
    print("Creating visualization...")
    print("="*70)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    batches = [r['batch_size'] for r in batch_results]
    accs = [r['test_acc'] for r in batch_results]
    ious = [r['test_iou'] for r in batch_results]
    maps = [r['test_map'] for r in batch_results]

    # Plot 1: Loss curves comparison
    colors = ['blue', 'orange', 'green', 'red']
    for idx, batch_size in enumerate(batch_sizes):
        if batch_size in training_histories:
            history = training_histories[batch_size]
            model = history['model']
            # Plot training losses
            if hasattr(model, 'train_losses') and len(model.train_losses) > 0:
                epochs = list(range(1, len(model.train_losses) + 1))
                axes[0, 0].plot(epochs, model.train_losses, 
                              label=f'Batch {batch_size}', color=colors[idx], linewidth=2)
    
    axes[0, 0].set_xlabel('Epoch', fontsize=12)
    axes[0, 0].set_ylabel('Training Loss', fontsize=12)
    axes[0, 0].set_title('Training Loss Curves', fontsize=14, fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Accuracy vs Batch Size
    axes[0, 1].plot(batches, accs, 'o-', linewidth=2, markersize=8)
    axes[0, 1].set_xlabel('Batch Size', fontsize=12)
    axes[0, 1].set_ylabel('Test Accuracy', fontsize=12)
    axes[0, 1].set_title('Accuracy vs Batch Size', fontsize=14, fontweight='bold')
    axes[0, 1].set_xticks(batches)
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: IoU vs Batch Size
    axes[1, 0].plot(batches, ious, 'o-', linewidth=2, markersize=8, color='orange')
    axes[1, 0].set_xlabel('Batch Size', fontsize=12)
    axes[1, 0].set_ylabel('Mean IoU', fontsize=12)
    axes[1, 0].set_title('IoU vs Batch Size', fontsize=14, fontweight='bold')
    axes[1, 0].set_xticks(batches)
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: mAP vs Batch Size
    axes[1, 1].plot(batches, maps, 'o-', linewidth=2, markersize=8, color='green')
    axes[1, 1].set_xlabel('Batch Size', fontsize=12)
    axes[1, 1].set_ylabel('mAP@0.5', fontsize=12)
    axes[1, 1].set_title('mAP vs Batch Size', fontsize=14, fontweight='bold')
    axes[1, 1].set_xticks(batches)
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("batch_size_results.png", dpi=150)
    plt.close()
    print("Visualization saved to batch_size_results.png")

    # Print summary
    print("\n" + "="*70)
    print("BATCH SIZE ANALYSIS SUMMARY")
    print("="*70)
    best_acc_idx = np.argmax(accs)
    best_iou_idx = np.argmax(ious)
    best_map_idx = np.argmax(maps)
    
    print("\nResults by batch size:")
    for r in batch_results:
        print(f"  Batch {r['batch_size']:3d}: Acc={r['test_acc']:.4f}, IoU={r['test_iou']:.4f}, mAP={r['test_map']:.4f}, Loss={r['test_loss']:.4f}")
    
    print(f"\nBest Accuracy: {accs[best_acc_idx]:.4f} with batch size {batches[best_acc_idx]}")
    print(f"Best IoU: {ious[best_iou_idx]:.4f} with batch size {batches[best_iou_idx]}")
    print(f"Best mAP: {maps[best_map_idx]:.4f} with batch size {batches[best_map_idx]}")
    print("\nNote: Larger batches provide stable gradients but may miss finer details.")
    print("="*70)
