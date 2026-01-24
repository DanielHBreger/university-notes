from torch.utils.data import DataLoader
import torch
import random
import numpy as np
import matplotlib.pyplot as plt

def train_one_epoch(model, dataloader, optimizer, device, cls_weight=1.0, box_weight=1.0, criterion_classification=None, criterion_bbox=None):
    model.train()
    running_total = 0.0
    running_cls = 0.0
    running_box = 0.0

    for images, labels, boxes, valids in dataloader:
        images = images.to(device)
        labels = labels.to(device)
        boxes = boxes.to(device)
        valids = valids.to(device)

        class_logits, bbox_pred = model(images)

        # classification loss (only on valid samples)
        if valids.sum() > 0:
            cls_loss = criterion_classification(class_logits[valids == 1], labels[valids == 1])
            box_loss = criterion_bbox(bbox_pred[valids == 1], boxes[valids == 1])
        else:
            cls_loss = torch.tensor(0.0, device=device)
            box_loss = torch.tensor(0.0, device=device)

        loss = cls_weight * cls_loss + box_weight * box_loss

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        running_total += loss.item()
        running_cls += cls_loss.item()
        running_box += box_loss.item()

    epoch_loss = running_total / max(1, len(dataloader))
    cls_loss = running_cls / max(1, len(dataloader))
    box_loss = running_box / max(1, len(dataloader))
    return epoch_loss, cls_loss, box_loss


def validate_one_epoch(model, dataloader, device, criterion_classification, criterion_bbox, cls_weight=1.0, box_weight=1.0,):
    model.eval()
    running_total = 0.0
    running_cls = 0.0
    running_box = 0.0

    with torch.no_grad():
        for images, labels, boxes, valids in dataloader:
            images = images.to(device)
            labels = labels.to(device)
            boxes = boxes.to(device)
            valids = valids.to(device)

            class_logits, bbox_pred = model(images)

            if valids.sum() > 0:
                cls_loss = criterion_classification(class_logits[valids == 1], labels[valids == 1])
                box_loss = criterion_bbox(bbox_pred[valids == 1], boxes[valids == 1])
            else:
                cls_loss = torch.tensor(0.0, device=device)
                box_loss = torch.tensor(0.0, device=device)

            loss = cls_weight * cls_loss + box_weight * box_loss

            running_total += loss.item()
            running_cls += cls_loss.item()
            running_box += box_loss.item()

    epoch_loss = running_total / max(1, len(dataloader))
    cls_loss = running_cls / max(1, len(dataloader))
    box_loss = running_box / max(1, len(dataloader))
    return epoch_loss, cls_loss, box_loss


def train_model(model, train_loader, valid_loader, optimizer, device, num_epochs, patience=10, criterion_classification=None, criterion_bbox=None):
    best_val_loss = float('inf')
    patience_counter = 0

    for epoch in range(num_epochs):
        train_loss, train_cls, train_box = train_one_epoch(
            model, train_loader, optimizer, device, cls_weight=1.0, box_weight=1.0, criterion_classification=criterion_classification, criterion_bbox=criterion_bbox
        )
        print(
            f"Epoch {epoch+1}/{num_epochs} | "
            f"Train: total={train_loss:.4f}, cls={train_cls:.4f}, box={train_box:.4f}"
        )

        val_loss, val_cls, val_box = validate_one_epoch(
            model, valid_loader, device, criterion_classification, criterion_bbox, cls_weight=1.0, box_weight=1.0
        )
        print(
            f"Epoch {epoch+1}/{num_epochs} | "
            f"Val: total={val_loss:.4f}, cls={val_cls:.4f}, box={val_box:.4f}"
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch+1}")
                break

def test_model(model, test_dataloader, device):
    # test the model
    test_loss, test_cls, test_box = validate_one_epoch(
        model, test_dataloader, device, cls_weight=1.0, box_weight=1.0
    )
    print(
        f"Test: total={test_loss:.4f}, cls={test_cls:.4f}, box={test_box:.4f}"
    )

def show_random(model, test_dataset, device, class_names, mean, std):
    # show a random sample prediction
    idx = random.randint(0, len(test_dataset) - 1)
    image, label, box, valid = test_dataset[idx]
    model.eval()
    with torch.no_grad():
        image = image.unsqueeze(0).to(device)
        class_logits, bbox_pred = model(image)
        predicted_label = torch.argmax(class_logits, dim=1).item()
        predicted_box = bbox_pred.squeeze(0).cpu().numpy()
    # visualize prediction
    image_np = image.squeeze(0).cpu().numpy().transpose(1, 2, 0)
    # denormalize
    image_np = (image_np * std) + mean
    image_np = np.clip(image_np, 0, 1)
    plt.imshow(image_np)
    ax = plt.gca()
    # ground truth box
    if valid.item() == 1:
        x_min, y_min, x_max, y_max = box.numpy()
        width = x_max - x_min
        height = y_max - y_min
        rect = plt.Rectangle((x_min, y_min), width, height, fill=False, color='green', linewidth=2)
        ax.add_patch(rect)
        ax.text(x_min, y_min, f"GT: {class_names[label.item()]}", fontsize=12, color='white', bbox=dict(facecolor='green', alpha=0.5))

    # predicted box
    x_min, y_min, x_max, y_max = predicted_box
    width = x_max - x_min
    height = y_max - y_min
    rect = plt.Rectangle((x_min, y_min), width, height, fill=False, color='red', linewidth=2)
    ax.add_patch(rect)
    ax.text(x_min, y_min + 15, f"Pred: {class_names[predicted_label]}", fontsize=12, color='white', bbox=dict(facecolor='red', alpha=0.5))
    plt.axis('off')
    plt.show()

def tranining_logic(train_dataset, valid_dataset, test_dataset, collate_fn, device, model, optimizer):
    # parameters
    batch_size = 800
    num_epochs = 20
    learning_rate = 1e-3
    num_workers = 4
    pin_memory = device == "cuda"
    persistent_workers = num_workers > 0

    train_dataloader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collate_fn,
        num_workers=num_workers,
        pin_memory=pin_memory
    )

    valid_dataloader = DataLoader(
        valid_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=num_workers,
        pin_memory=pin_memory
    )

    test_dataloader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=num_workers,
        pin_memory=pin_memory
    )
    train_model(model, train_dataloader, valid_dataloader, optimizer, device, num_epochs=num_epochs)